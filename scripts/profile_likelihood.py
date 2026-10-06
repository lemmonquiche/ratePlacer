#!/usr/bin/env python3
"""
Profile log-likelihood estimation of alpha with confidence intervals.

The observed D^(k) matrices are treated as empirical frequencies from
multinomial sampling.  Counts are derived from a user-specified number
of reads and read length.
"""

import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.linalg import inv
import matplotlib.pyplot as plt
from dataclasses import dataclass
from damage_estimator import (
    compute_A, compute_alpha_max, compute_B_analytic, find_feasible_interval,
    _make_row_sum_constraints, _initial_guess
)


@dataclass
class ProfileLikelihoodResult:
    alpha_hat: float
    B_hat: np.ndarray
    loglik_max: float
    ci_lo: float
    ci_hi: float
    feasible_interval: tuple
    profile_alphas: np.ndarray
    profile_logliks: np.ndarray
    counts: np.ndarray  # (K, 4, 4)
    n_per_category: np.ndarray  # (K,) total counts per category
    alpha_upper: float = None  # upper end of the alpha search
    ci_hi_censored: bool = False  # True if the profile never fell below the CI threshold


# Default nucleotide frequencies (A, C, G, T)
DEFAULT_BASE_FREQ = np.array([0.32461449, 0.17526887, 0.17534515, 0.32477149])


def resolve_alpha_upper(Q, alpha_upper=None, verbose=True):
    """Upper end of the alpha search.

    Defaults to alpha_max = 1/max_i(-Q_ii), the largest valid alpha (Section 4.3.1).
    A user-supplied alpha_upper can lower this bound but not raise it.
    """
    alpha_max = compute_alpha_max(Q)
    if alpha_upper is None:
        return alpha_max
    if alpha_upper > alpha_max:
        if verbose:
            print(f"  alpha_upper={alpha_upper:g} exceeds alpha_max={alpha_max:.6g}; using alpha_max")
        return alpha_max
    return alpha_upper


def alpha_grid(alpha_upper, n_grid):
    """Profile scan grid on [0, alpha_upper]: 0 plus n_grid-1 log-spaced points.

    Log spacing resolves the small alphas typical of ancient samples while still
    covering the full range up to alpha_max.
    """
    return np.concatenate([[0.0], np.geomspace(alpha_upper * 1e-6, alpha_upper, n_grid - 1)])


def _profile_ci(profile_ll, alpha_upper, n_grid, ci_level, verbose):
    """Grid scan, MLE refinement and likelihood-ratio CI for a profile log-likelihood.

    profile_ll(alpha, B_prev) returns (loglik, B).  The CI bounds are found by walking
    outward from the MLE to the first grid point below the threshold and root-finding
    in between.  If the profile never falls below the threshold above the MLE, the
    upper bound is censored at alpha_upper.

    Returns (alpha_hat, ll_max, B_hat, ci_lo, ci_hi, ci_hi_censored, alphas, logliks).
    """
    alphas = alpha_grid(alpha_upper, n_grid)
    logliks = np.empty(len(alphas))
    B_carry = None

    if verbose:
        print(f"Profile scan over [0, {alpha_upper:.6g}] with {len(alphas)} points ...")

    for i, alpha in enumerate(alphas):
        logliks[i], B_carry = profile_ll(alpha, B_carry)

    best_idx = int(np.argmax(logliks))
    if verbose:
        print(f"  Grid max at alpha = {alphas[best_idx]:.6g}, loglik = {logliks[best_idx]:.4f}")

    # Refine alpha_hat between the neighbouring grid points
    bracket_lo = alphas[max(best_idx - 2, 0)]
    bracket_hi = alphas[min(best_idx + 2, len(alphas) - 1)]
    ref = minimize_scalar(lambda a: -profile_ll(a, None)[0], bounds=(bracket_lo, bracket_hi),
                          method='bounded', options={'xatol': 1e-10, 'maxiter': 200})
    alpha_hat, ll_max = ref.x, -ref.fun
    if logliks[best_idx] > ll_max:
        alpha_hat, ll_max = alphas[best_idx], logliks[best_idx]

    _, B_hat = profile_ll(alpha_hat, None)

    if verbose:
        print(f"  Refined: alpha = {alpha_hat:.8f}, loglik = {ll_max:.4f}")

    threshold = ll_max - ci_level

    def crossing(a_in, a_out):
        """Alpha between a_in (above the threshold) and a_out (below it) where the profile crosses it."""
        try:
            return brentq(lambda a: profile_ll(a, None)[0] - threshold,
                          min(a_in, a_out), max(a_in, a_out), xtol=1e-12)
        except ValueError:
            if verbose:
                print(f"  Note: no sign change between {a_in:.6g} and {a_out:.6g}; using {a_out:.6g}")
            return a_out

    # Upper bound: first grid point above alpha_hat that is below the threshold
    ci_hi, ci_hi_censored = alpha_upper, True
    prev = alpha_hat
    for a, ll in zip(alphas, logliks):
        if a <= alpha_hat:
            continue
        if ll < threshold:
            ci_hi, ci_hi_censored = crossing(prev, a), False
            break
        prev = a

    # Lower bound: first grid point below alpha_hat that is below the threshold (else 0)
    ci_lo = 0.0
    prev = alpha_hat
    for a, ll in zip(alphas[::-1], logliks[::-1]):
        if a >= alpha_hat:
            continue
        if ll < threshold:
            ci_lo = crossing(prev, a)
            break
        prev = a

    if verbose:
        print(f"  95% CI: [{ci_lo:.8f}, {ci_hi:.8f}]")
        if ci_hi_censored:
            print(f"  WARNING: the profile log-likelihood stays within {ci_level} of its maximum "
                  f"up to alpha = {alpha_upper:.6g}, so the upper CI bound is censored there.")

    return alpha_hat, ll_max, B_hat, ci_lo, ci_hi, ci_hi_censored, alphas, logliks


def compute_counts(D, n_reads, read_length, base_freq=None):
    """Derive expected counts from D matrices, number of reads, and read length.

    Categories: positions 1..n_end from 5', middle, positions 1..n_end from 3'.
    Each endpoint position gets n_reads observations.
    Middle gets n_reads * (read_length - 2*n_end) observations.

    Parameters
    ----------
    base_freq : (4,) nucleotide frequencies for A, C, G, T.
        n^(k)_i = n^(k) * base_freq[i].  If None, uses DEFAULT_BASE_FREQ.

    Returns
    -------
    counts : (K, 4, 4) expected count matrix
    n_per_cat : (K,) total observations per category
    """
    if base_freq is None:
        base_freq = DEFAULT_BASE_FREQ
    base_freq = np.asarray(base_freq)

    K = D.shape[0]
    n_end = (K - 1) // 2
    n_middle_positions = read_length - 2 * n_end

    n_per_cat = np.empty(K)
    n_per_cat[:n_end] = n_reads
    n_per_cat[n_end] = n_reads * n_middle_positions
    n_per_cat[n_end + 1:] = n_reads

    # Counts: c^(k)_{ij} = n^(k) * base_freq[i] * D^(k)_{ij}
    counts = np.zeros_like(D)
    for k in range(K):
        for i in range(4):
            counts[k, i, :] = n_per_cat[k] * base_freq[i] * D[k, i, :]

    return counts, n_per_cat


def loglik_multinomial(D, C, B, A, n_per_cat, base_freq=None):
    """Compute the multinomial log-likelihood.

    L = sum_k sum_i n^(k) * base_freq[i] * sum_j D^(k)_{ij} log(P^(k)_{ij})
    where P^(k) = A B^(k) C.
    """
    if base_freq is None:
        base_freq = DEFAULT_BASE_FREQ
    K = D.shape[0]
    eps = 1e-300
    ll = 0.0
    for k in range(K):
        P = A @ B[k] @ C
        P = np.maximum(P, eps)
        for i in range(4):
            ll += n_per_cat[k] * base_freq[i] * np.sum(D[k, i, :] * np.log(P[i, :]))
    return ll


def solve_inner_mll_counts(D_k, C, A, n_k, B0=None, base_freq=None):
    """Solve max sum_i n_k*f_i sum_j D_k_ij log(A B C)_ij  s.t. B >= 0, rows sum to 1.

    Uses per-nucleotide frequency weights f_i.  The optimum B does not depend
    on the overall scale n_k, only on f_i-weighted cross-entropy, so we optimise
    with weights f_i for numerical stability and scale by n_k at the end.

    Returns (B, loglik_contribution).
    """
    if base_freq is None:
        base_freq = DEFAULT_BASE_FREQ

    M = np.kron(A, C.T)
    d = D_k.ravel()  # (16,)

    # Build per-element weight vector: w[4*i+j] = base_freq[i]
    w = np.repeat(base_freq, 4)  # (16,)

    eps = 1e-300

    # Optimize with f_i weights (no n_k scaling)
    def neg_loglik_unit(b):
        p = M @ b
        p = np.maximum(p, eps)
        return -(w * d) @ np.log(p)

    def neg_loglik_unit_grad(b):
        p = M @ b
        p = np.maximum(p, eps)
        return -M.T @ (w * d / p)

    constraints = _make_row_sum_constraints()
    bounds = [(0.0, None)] * 16

    # Try multiple starting points, keep best
    candidates = []
    if B0 is not None:
        candidates.append(B0)
    candidates.append(_initial_guess(D_k, C, A))
    candidates.append(np.full((4, 4), 0.25))  # uniform

    best_b = None
    best_val = np.inf

    for B_init in candidates:
        b0 = B_init.ravel()
        res = minimize(neg_loglik_unit, b0, jac=neg_loglik_unit_grad, method='SLSQP',
                       bounds=bounds, constraints=constraints,
                       options={'maxiter': 5000, 'ftol': 1e-15})
        if res.fun < best_val:
            best_val = res.fun
            best_b = res.x

    B = best_b.reshape(4, 4)
    ll = -n_k * best_val  # scale by n_k
    return B, ll


def profile_loglik(alpha, Q, C, D, n_per_cat, B_prev=None, base_freq=None):
    """Compute the profile log-likelihood at a given alpha."""
    A = compute_A(alpha, Q)
    K = D.shape[0]
    B_all = np.empty((K, 4, 4))
    total_ll = 0.0

    for k in range(K):
        B0 = B_prev[k] if B_prev is not None else None
        B_k, ll_k = solve_inner_mll_counts(D[k], C, A, n_per_cat[k], B0=B0,
                                           base_freq=base_freq)
        B_all[k] = B_k
        total_ll += ll_k

    return total_ll, B_all


def estimate_alpha_mll(Q, C, D, n_reads=1000, read_length=50,
                       alpha_upper=None, n_grid=100, ci_level=1.92,
                       base_freq=None, verbose=True):
    """Profile likelihood estimation of alpha with CI.

    Parameters
    ----------
    Q : (4,4) rate matrix
    C : (4,4) sequencing error matrix
    D : (K,4,4) observed mismatch frequencies
    n_reads : number of reads
    read_length : read length in bp
    alpha_upper : upper bound for search (default and maximum: alpha_max = 1/max_i(-Q_ii))
    n_grid : grid points for profile scan
    ci_level : log-likelihood drop for CI (1.92 for 95% CI)
    base_freq : (4,) nucleotide frequencies [A, C, G, T], or None for default
    verbose : print progress

    Returns
    -------
    ProfileLikelihoodResult
    """
    if base_freq is None:
        base_freq = DEFAULT_BASE_FREQ
    K = D.shape[0]
    counts, n_per_cat = compute_counts(D, n_reads, read_length, base_freq=base_freq)
    alpha_upper = resolve_alpha_upper(Q, alpha_upper, verbose)

    if verbose:
        print(f"Counts: {n_reads} reads, {read_length} bp")
        print(f"  Endpoint categories: {n_per_cat[0]:.0f} obs each")
        print(f"  Middle category: {n_per_cat[K//2]:.0f} obs")

    # Find feasible interval for reference
    feas = find_feasible_interval(Q, C, D, alpha_upper=alpha_upper)
    feas_lo = feas[0] if feas[0] is not None else 0.0
    feas_hi = feas[1] if feas[1] is not None else 0.0
    if verbose:
        print(f"  Feasible interval: [{feas_lo:.8f}, {feas_hi:.8f}]")

    def profile_ll(alpha, B_prev):
        return profile_loglik(alpha, Q, C, D, n_per_cat, B_prev, base_freq=base_freq)

    alpha_hat, ll_max, B_hat, ci_lo, ci_hi, ci_hi_censored, alphas, logliks = _profile_ci(
        profile_ll, alpha_upper, n_grid, ci_level, verbose)

    return ProfileLikelihoodResult(
        alpha_hat=alpha_hat,
        B_hat=B_hat,
        loglik_max=ll_max,
        ci_lo=ci_lo,
        ci_hi=ci_hi,
        feasible_interval=(feas_lo, feas_hi),
        profile_alphas=alphas,
        profile_logliks=logliks,
        counts=counts,
        n_per_category=n_per_cat,
        alpha_upper=alpha_upper,
        ci_hi_censored=ci_hi_censored,
    )


def plot_profile_likelihood(result, title='', save_path=None):
    """Plot profile log-likelihood with CI."""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(result.profile_alphas, result.profile_logliks, 'b-', lw=1.5)
    ax.axvline(result.alpha_hat, color='red', ls='--', lw=1,
               label=f'MLE: {result.alpha_hat:.6f}')

    # CI
    threshold = result.loglik_max - 1.92
    ax.axhline(threshold, color='gray', ls=':', lw=0.8, label='95% CI threshold')
    ax.axvspan(result.ci_lo, result.ci_hi, alpha=0.12, color='red',
               label=f'95% CI: [{result.ci_lo:.4f}, {result.ci_hi:.4f}]')

    # Feasible interval
    flo, fhi = result.feasible_interval
    if fhi > flo:
        ax.axvspan(flo, fhi, alpha=0.12, color='green', label=f'Feasible interval')

    # Zoom to the CI: the scan runs to alpha_max, far beyond the alphas of interest
    alphas = np.asarray(result.profile_alphas)
    logliks = np.asarray(result.profile_logliks)
    span = max(result.ci_hi - result.ci_lo, result.ci_hi, 1e-12)
    x_hi = min(alphas[-1], result.ci_hi + span)
    visible = alphas <= x_hi
    y_lo = min(logliks[visible].min(), threshold)
    y_pad = 0.05 * (result.loglik_max - y_lo) or 1.0
    ax.set_xlim(0, x_hi)
    ax.set_ylim(y_lo - y_pad, result.loglik_max + y_pad)

    ax.set_xlabel('alpha')
    ax.set_ylabel('Profile log-likelihood')
    ax.set_title(title or 'Profile log-likelihood')
    ax.legend(fontsize=9)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"  Saved {save_path}")
    return fig, ax


def solve_inner_mll_counts_ref(D_k, C, A, n_ref_k, B0=None):
    """Solve max sum_{i,j} n_ref_k[i] * D_k[i,j] * log(A B C)_{i,j}  s.t. B >= 0, rows sum to 1.

    Parameters
    ----------
    n_ref_k : (4,) array of reference base counts [A, C, G, T] for this category.
        Encodes both total coverage (n_ref_k.sum()) and per-nucleotide composition.

    Returns (B, loglik_contribution).
    """
    M = np.kron(A, C.T)
    d = D_k.ravel()

    n_k = n_ref_k.sum()
    base_freq_k = n_ref_k / n_k  # per-category empirical base frequencies
    w = np.repeat(base_freq_k, 4)  # w[4*i+j] = base_freq_k[i]

    eps = 1e-300

    def neg_loglik_unit(b):
        p = M @ b
        p = np.maximum(p, eps)
        return -(w * d) @ np.log(p)

    def neg_loglik_unit_grad(b):
        p = M @ b
        p = np.maximum(p, eps)
        return -M.T @ (w * d / p)

    constraints = _make_row_sum_constraints()
    bounds = [(0.0, None)] * 16

    candidates = []
    if B0 is not None:
        candidates.append(B0)
    candidates.append(_initial_guess(D_k, C, A))
    candidates.append(np.full((4, 4), 0.25))

    best_b = None
    best_val = np.inf

    for B_init in candidates:
        b0 = B_init.ravel()
        res = minimize(neg_loglik_unit, b0, jac=neg_loglik_unit_grad, method='SLSQP',
                       bounds=bounds, constraints=constraints,
                       options={'maxiter': 5000, 'ftol': 1e-15})
        if res.fun < best_val:
            best_val = res.fun
            best_b = res.x

    B = best_b.reshape(4, 4)
    ll = -n_k * best_val
    return B, ll


def profile_loglik_ref(alpha, Q, C, D, n_ref, B_prev=None):
    """Profile log-likelihood at alpha using per-category reference counts.

    Parameters
    ----------
    n_ref : (K, 4) array; n_ref[k, i] = reference count of nucleotide i at category k.
    """
    A = compute_A(alpha, Q)
    K = D.shape[0]
    B_all = np.empty((K, 4, 4))
    total_ll = 0.0

    for k in range(K):
        B0 = B_prev[k] if B_prev is not None else None
        B_k, ll_k = solve_inner_mll_counts_ref(D[k], C, A, n_ref[k], B0=B0)
        B_all[k] = B_k
        total_ll += ll_k

    return total_ll, B_all


def estimate_alpha_mll_ref(Q, C, D, n_ref, alpha_upper=None, n_grid=100,
                            ci_level=1.92, verbose=True):
    """Profile likelihood estimation of alpha using per-category reference counts.

    Drop-in alternative to estimate_alpha_mll that accepts actual reference counts
    instead of (n_reads, read_length, base_freq).  This avoids two approximations
    in the original: (1) equal depth across endpoint positions, and (2) global
    base frequencies applied uniformly to all position categories.

    Parameters
    ----------
    Q : (4, 4) rate matrix
    C : (4, 4) sequencing error matrix
    D : (K, 4, 4) observed mismatch frequency matrices
    n_ref : (K, 4) reference count array; n_ref[k, i] = number of reference bases
        of type i observed at position category k.  Encodes both per-category
        coverage and per-category nucleotide composition.
    alpha_upper : upper bound for alpha search and feasibility scan
        (default and maximum: alpha_max = 1/max_i(-Q_ii))
    n_grid : number of grid points for the coarse profile scan
    ci_level : log-likelihood drop for CI (default 1.92 gives 95% CI)
    verbose : print progress

    Returns
    -------
    ProfileLikelihoodResult
        Same dataclass as estimate_alpha_mll.  The counts field holds the implied
        count matrix n_ref[k, i] * D[k, i, j]; n_per_category holds n_ref.sum(axis=1).
    """
    n_ref = np.asarray(n_ref, dtype=float)
    K = D.shape[0]
    if n_ref.shape != (K, 4):
        raise ValueError(f"n_ref must have shape (K, 4) = ({K}, 4), got {n_ref.shape}")

    n_per_cat = n_ref.sum(axis=1)
    alpha_upper = resolve_alpha_upper(Q, alpha_upper, verbose)

    if verbose:
        print(f"Reference counts: K={K} categories")
        print(f"  Per-category total observations: min={n_per_cat.min():.0f}, "
              f"max={n_per_cat.max():.0f}, sum={n_per_cat.sum():.0f}")

    feas = find_feasible_interval(Q, C, D, alpha_upper=alpha_upper)
    feas_lo = feas[0] if feas[0] is not None else 0.0
    feas_hi = feas[1] if feas[1] is not None else 0.0
    if verbose:
        print(f"  Feasible interval: [{feas_lo:.8f}, {feas_hi:.8f}]")

    def profile_ll(alpha, B_prev):
        return profile_loglik_ref(alpha, Q, C, D, n_ref, B_prev)

    alpha_hat, ll_max, B_hat, ci_lo, ci_hi, ci_hi_censored, alphas, logliks = _profile_ci(
        profile_ll, alpha_upper, n_grid, ci_level, verbose)

    counts = np.zeros((K, 4, 4))
    for k in range(K):
        for i in range(4):
            counts[k, i, :] = n_ref[k, i] * D[k, i, :]

    return ProfileLikelihoodResult(
        alpha_hat=alpha_hat,
        B_hat=B_hat,
        loglik_max=ll_max,
        ci_lo=ci_lo,
        ci_hi=ci_hi,
        feasible_interval=(feas_lo, feas_hi),
        profile_alphas=alphas,
        profile_logliks=logliks,
        counts=counts,
        n_per_category=n_per_cat,
        alpha_upper=alpha_upper,
        ci_hi_censored=ci_hi_censored,
    )


if __name__ == '__main__':
    import os
    DATADIR = os.path.dirname(os.path.abspath(__file__))
    OUTDIR = os.path.join(DATADIR, 'data_analysis_results')

    for prefix in ['ss', 'ds']:
        Q = np.load(os.path.join(DATADIR, f'{prefix}_A.npy'))
        C = np.load(os.path.join(DATADIR, f'{prefix}_C.npy'))
        D = np.load(os.path.join(DATADIR, f'{prefix}_D.npy'))

        print(f"\n{'#'*60}")
        print(f"# {prefix.upper()} — Profile likelihood")
        print(f"{'#'*60}")

        result = estimate_alpha_mll(Q, C, D, n_reads=1000, read_length=50,
                                    n_grid=100, verbose=True)

        plot_profile_likelihood(
            result,
            title=f'{prefix.upper()} — Profile log-likelihood (n=1000, L=50bp)',
            save_path=os.path.join(OUTDIR, f'{prefix}_profile_likelihood.png'))

        print(f"\n  Summary:")
        print(f"    alpha_hat = {result.alpha_hat:.8f}")
        print(f"    95% CI = [{result.ci_lo:.8f}, {result.ci_hi:.8f}]")
        print(f"    Feasible interval = [{result.feasible_interval[0]:.8f}, "
              f"{result.feasible_interval[1]:.8f}]")
