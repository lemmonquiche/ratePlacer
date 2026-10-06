"""
Damage matrix estimation via constrained matrix factorisation.

Implements the model D^(k) = A(alpha) B^(k) C where:
  D^(k) : observed mismatch matrices (K x 4 x 4), known
  C      : sequencing error matrix (4 x 4), known
  A(a)   : expm(alpha * Q), evolutionary substitution (4 x 4), known up to scalar alpha
  B^(k)  : damage matrices (K x 4 x 4), estimated

All matrices are row-stochastic, M[from, to], so the product reads left to right in the
physical order reference -> evolution (A) -> damage (B) -> sequencing error (C) -> read.

Two loss functions: weighted least squares (WLS) and multinomial log-likelihood (MLL).
Estimation via profile optimisation over alpha with constrained inner solves for B^(k).
"""

import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq
from scipy.linalg import inv, expm
import matplotlib.pyplot as plt
from dataclasses import dataclass


@dataclass
class EstimationResult:
    alpha_hat: float
    B_hat: np.ndarray
    loss: float
    feasible_interval: tuple
    profile_alphas: np.ndarray
    profile_losses: np.ndarray
    reconstruction_errors: np.ndarray
    loss_type: str


def default_weights(read_length=50, K=21):
    """Default category weights assuming K = 2*n_end + 1 categories.

    The first n_end categories are 5' end positions, then 1 middle category,
    then n_end categories for 3' end positions.  Each endpoint position gets
    weight 1/read_length; the middle gets the remainder.
    """
    n_end = (K - 1) // 2
    pi = np.empty(K)
    w_end = 1.0 / read_length
    pi[:n_end] = w_end
    pi[n_end] = 1.0 - 2 * n_end * w_end
    pi[n_end + 1:] = w_end
    return pi


def compute_A(alpha, Q):
    """Substitution matrix A(alpha) = expm(alpha * Q)."""
    return expm(alpha * Q)


def compute_alpha_max(Q):
    """Largest valid alpha, 1 / max_i(-Q_ii) (Section 4.3.1 of the paper)."""
    return 1.0 / np.max(-np.diag(Q))


def compute_B_analytic(alpha, Q, C, D):
    """Analytic (unconstrained) B^(k) = A_inv @ D^(k) @ C_inv for all k."""
    C_inv = inv(C)
    A_inv = inv(compute_A(alpha, Q))
    K = D.shape[0]
    B = np.empty_like(D)
    for k in range(K):
        B[k] = A_inv @ D[k] @ C_inv
    return B


def find_feasible_interval(Q, C, D, alpha_upper=2.0, n_grid=2000, tol=0.0):
    """Find the interval of alpha where all analytic B^(k) are non-negative.

    Returns (alpha_lo, alpha_hi).  If no feasible point exists, returns (None, None).
    """
    alphas = np.linspace(0, alpha_upper, n_grid)

    def min_B_entry(alpha):
        B = compute_B_analytic(alpha, Q, C, D)
        return B.min()

    # Evaluate on grid
    min_vals = np.array([min_B_entry(a) for a in alphas])
    feasible = min_vals >= tol

    if not feasible.any():
        return (None, None)

    idx = np.where(feasible)[0]
    lo_idx, hi_idx = idx[0], idx[-1]

    # Refine boundaries with root-finding
    alpha_lo = alphas[lo_idx]
    alpha_hi = alphas[hi_idx]

    if lo_idx > 0:
        try:
            alpha_lo = brentq(lambda a: min_B_entry(a) - tol,
                              alphas[lo_idx - 1], alphas[lo_idx], xtol=1e-12)
        except ValueError:
            pass

    if hi_idx < n_grid - 1:
        try:
            alpha_hi = brentq(lambda a: min_B_entry(a) - tol,
                              alphas[hi_idx], alphas[hi_idx + 1], xtol=1e-12)
        except ValueError:
            pass

    return (alpha_lo, alpha_hi)


def _make_row_sum_constraints(n=4):
    """Equality constraints: each row of a 4x4 matrix (flattened to 16-vec) sums to 1."""
    constraints = []
    for i in range(n):
        idx = list(range(i * n, (i + 1) * n))
        constraints.append({
            'type': 'eq',
            'fun': lambda b, idx=idx: np.sum(b[idx]) - 1.0,
            'jac': lambda b, idx=idx: _row_sum_jac(idx, n * n)
        })
    return constraints


def _row_sum_jac(idx, size):
    j = np.zeros(size)
    j[idx] = 1.0
    return j


def _initial_guess(D_k, C, A):
    """Warm-start: analytic inverse clamped to non-negative and row-normalised."""
    try:
        B0 = inv(A) @ D_k @ inv(C)
    except np.linalg.LinAlgError:
        B0 = np.full((4, 4), 0.25)
    B0 = np.maximum(B0, 1e-10)
    B0 /= B0.sum(axis=1, keepdims=True)
    return B0


def solve_inner_wls(D_k, C, A, B0=None):
    """Solve min ||A @ B @ C - D_k||_F^2  s.t. B >= 0, rows sum to 1.

    Uses Kronecker vectorisation: for row-major vec, vec(A B C) = kron(A, C^T) vec(B).
    """
    M = np.kron(A, C.T)  # (16, 16)
    d = D_k.ravel()       # (16,)

    if B0 is None:
        B0 = _initial_guess(D_k, C, A)
    b0 = B0.ravel()

    def obj(b):
        r = M @ b - d
        return 0.5 * r @ r

    def grad(b):
        return M.T @ (M @ b - d)

    constraints = _make_row_sum_constraints()
    bounds = [(0.0, None)] * 16

    res = minimize(obj, b0, jac=grad, method='SLSQP',
                   bounds=bounds, constraints=constraints,
                   options={'maxiter': 500, 'ftol': 1e-15})
    B = res.x.reshape(4, 4)
    loss = 2.0 * res.fun  # ||.||_F^2
    return B, loss


def solve_inner_mll(D_k, C, A, B0=None):
    """Solve min -sum D_k_ij * log(A B C)_ij  s.t. B >= 0, rows sum to 1."""
    M = np.kron(A, C.T)
    d = D_k.ravel()

    if B0 is None:
        B0 = _initial_guess(D_k, C, A)
    b0 = B0.ravel()

    eps = 1e-15

    def obj(b):
        p = M @ b
        p = np.maximum(p, eps)
        return -d @ np.log(p)

    def grad(b):
        p = M @ b
        p = np.maximum(p, eps)
        return -M.T @ (d / p)

    constraints = _make_row_sum_constraints()
    bounds = [(0.0, None)] * 16

    res = minimize(obj, b0, jac=grad, method='SLSQP',
                   bounds=bounds, constraints=constraints,
                   options={'maxiter': 1000, 'ftol': 1e-15})
    B = res.x.reshape(4, 4)
    loss = res.fun
    return B, loss


def solve_all_inner(D, C, A, pi, loss_type='WLS', B_prev=None):
    """Solve all K inner problems.  Returns (B_all, total_weighted_loss)."""
    K = D.shape[0]
    B_all = np.empty((K, 4, 4))
    total_loss = 0.0
    solver = solve_inner_wls if loss_type == 'WLS' else solve_inner_mll

    for k in range(K):
        B0 = B_prev[k] if B_prev is not None else None
        B_k, loss_k = solver(D[k], C, A, B0=B0)
        B_all[k] = B_k
        total_loss += pi[k] * loss_k

    return B_all, total_loss


def profile_objective(alpha, Q, C, D, pi, loss_type='WLS', B_prev=None):
    """Compute L*(alpha) = min_B L(alpha, B)."""
    A = compute_A(alpha, Q)
    B_all, loss = solve_all_inner(D, C, A, pi, loss_type, B_prev)
    return loss, B_all


def estimate_damage(Q, C, D, loss_type='WLS', pi=None, n_grid=50,
                    alpha_range=None, alpha_upper=2.0, verbose=True):
    """Main entry point: profile optimisation over alpha.

    Parameters
    ----------
    Q : (4, 4) generator matrix
    C : (4, 4) sequencing error matrix
    D : (K, 4, 4) observed mismatch matrices
    loss_type : 'WLS' or 'MLL'
    pi : (K,) category weights, or None for default (50bp read assumption)
    n_grid : number of grid points for the coarse search
    alpha_range : (lo, hi) to override the search range, or None for auto
    alpha_upper : upper bound for feasibility scan and default search range
    verbose : print progress

    Returns
    -------
    EstimationResult
    """
    K = D.shape[0]
    if pi is None:
        pi = default_weights(read_length=50, K=K)
    pi = np.asarray(pi, dtype=float)
    pi /= pi.sum()

    # Step 1: find feasible interval
    if verbose:
        print(f"Finding feasible interval (alpha_upper = {alpha_upper:.6f}) ...")
    feas = find_feasible_interval(Q, C, D, alpha_upper=alpha_upper)
    if verbose:
        if feas[0] is not None:
            print(f"  Feasible interval: [{feas[0]:.8f}, {feas[1]:.8f}]  "
                  f"(width = {feas[1] - feas[0]:.2e})")
        else:
            print("  No exact feasible interval found (will use full range with constrained solver)")

    # Determine search range
    if alpha_range is not None:
        lo, hi = alpha_range
    elif feas[0] is not None:
        margin = 0.05 * (feas[1] - feas[0] + 1e-8)
        lo = max(0, feas[0] - margin)
        hi = feas[1] + margin
    else:
        lo, hi = 0.0, alpha_upper

    # Step 2: grid search
    if verbose:
        print(f"Grid search over alpha in [{lo:.6f}, {hi:.6f}] with {n_grid} points ...")
    alphas = np.linspace(lo, hi, n_grid)
    losses = np.empty(n_grid)
    B_carry = None

    for i, alpha in enumerate(alphas):
        loss_val, B_carry = profile_objective(alpha, Q, C, D, pi, loss_type, B_carry)
        losses[i] = loss_val

    best_idx = np.argmin(losses)
    if verbose:
        print(f"  Grid minimum at alpha = {alphas[best_idx]:.8f}, loss = {losses[best_idx]:.8e}")

    # Step 3: refine with bounded scalar minimisation
    bracket_lo = alphas[max(best_idx - 1, 0)]
    bracket_hi = alphas[min(best_idx + 1, n_grid - 1)]

    if verbose:
        print(f"Refining alpha in [{bracket_lo:.8f}, {bracket_hi:.8f}] ...")

    def scalar_obj(alpha):
        loss_val, _ = profile_objective(alpha, Q, C, D, pi, loss_type)
        return loss_val

    ref = minimize_scalar(scalar_obj, bounds=(bracket_lo, bracket_hi), method='bounded',
                          options={'xatol': 1e-12, 'maxiter': 200})
    alpha_hat = ref.x

    # Step 4: final solve
    _, B_hat = profile_objective(alpha_hat, Q, C, D, pi, loss_type)
    final_loss = ref.fun

    if verbose:
        print(f"  Refined alpha = {alpha_hat:.10f}, loss = {final_loss:.8e}")

    # Reconstruction errors
    A_hat = compute_A(alpha_hat, Q)
    recon_errors = np.array([
        np.linalg.norm(A_hat @ B_hat[k] @ C - D[k], 'fro') for k in range(K)
    ])

    if verbose:
        print(f"  Mean reconstruction error (Frobenius): {recon_errors.mean():.6e}")
        print(f"  Max  reconstruction error (Frobenius): {recon_errors.max():.6e}")

    return EstimationResult(
        alpha_hat=alpha_hat,
        B_hat=B_hat,
        loss=final_loss,
        feasible_interval=feas,
        profile_alphas=alphas,
        profile_losses=losses,
        reconstruction_errors=recon_errors,
        loss_type=loss_type,
    )


def plot_profile(result, title='', save_path=None):
    """Plot the profile objective L*(alpha) vs alpha."""
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(result.profile_alphas, result.profile_losses, 'b-', lw=1.5)
    ax.axvline(result.alpha_hat, color='r', ls='--', lw=1, label=f'alpha = {result.alpha_hat:.6f}')

    if result.feasible_interval[0] is not None:
        ax.axvspan(result.feasible_interval[0], result.feasible_interval[1],
                   alpha=0.15, color='green', label='Feasible interval')

    ax.set_xlabel('alpha')
    ax.set_ylabel(f'Profile loss ({result.loss_type})')
    ax.set_title(title or f'Profile optimisation ({result.loss_type})')
    ax.legend()
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"  Profile plot saved to {save_path}")
    return fig, ax


def print_diagnostics(result, nucleotides='ACGT'):
    """Print summary diagnostics."""
    print(f"\n{'='*60}")
    print(f"Loss type: {result.loss_type}")
    print(f"Estimated alpha: {result.alpha_hat:.10f}")
    if result.feasible_interval[0] is not None:
        lo, hi = result.feasible_interval
        print(f"Feasible interval: [{lo:.8f}, {hi:.8f}]  (width = {hi - lo:.2e})")
    else:
        print("Feasible interval: none (constrained solver used)")
    print(f"Final loss: {result.loss:.8e}")
    print(f"Reconstruction error (mean): {result.reconstruction_errors.mean():.6e}")
    print(f"Reconstruction error (max):  {result.reconstruction_errors.max():.6e}")

    K = result.B_hat.shape[0]
    n_end = (K - 1) // 2

    # Print damage matrix for first 5' position, middle, and first 3' position
    for label, idx in [("5' pos 1", 0), ("middle", n_end), ("3' pos 1", n_end + 1)]:
        B = result.B_hat[idx]
        print(f"\nB^({label}) [category {idx}]:")
        header = "     " + "  ".join(f"{n:>8s}" for n in nucleotides)
        print(header)
        for i, ni in enumerate(nucleotides):
            vals = "  ".join(f"{B[i, j]:8.5f}" for j in range(4))
            print(f"  {ni}  {vals}")

    print(f"{'='*60}\n")
