#!/usr/bin/env python3
"""
Run full damage estimation analysis for ss and ds data:
  1. Feasible interval for alpha (searched up to alpha_max = 1/max_i(-Q_ii))
  2. Profile likelihood with confidence interval
  3. B matrices at the CI boundaries
  4. Damage profile plots
  5. All results saved to data_analysis_results/
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import inv
from damage_estimator import (
    compute_A, compute_alpha_max, compute_B_analytic, find_feasible_interval, default_weights
)
from profile_likelihood import (
    estimate_alpha_mll, estimate_alpha_mll_ref,
    plot_profile_likelihood,
    profile_loglik, profile_loglik_ref,
    compute_counts,
)

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_analysis_results')
DATADIR = os.path.dirname(os.path.abspath(__file__))

# User parameters (per dataset)
PARAMS = {
    'ss': {'n_reads': 95202, 'read_length': 36.4},
    'ds': {'n_reads': 152192, 'read_length': 35.6},
}


def total_damage_rate(B):
    """Per-category average off-diagonal probability (total damage/error rate)."""
    K = B.shape[0]
    rates = np.zeros(K)
    for k in range(K):
        rates[k] = 1.0 - np.mean(np.diag(B[k]))
    return rates


def ct_damage_rate(B):
    """Per-category C->T damage rate (row 1 = C, col 3 = T in ACGT order)."""
    return B[:, 1, 3]


def position_coords(K=21, read_length=50):
    """X-coordinates for plotting: 5' positions, gap, middle, gap, 3' positions."""
    n_end = (K - 1) // 2
    rl = int(round(read_length))
    x = list(range(1, n_end + 1))
    x.append(rl // 2)
    x += list(range(rl - n_end + 1, rl + 1))
    return np.array(x)


def plot_damage_profiles(B_lo, B_hi, prefix, alpha_lo, alpha_hi, read_length, out_dir):
    """Plot total damage rate and C->T rate along sequence positions."""
    K = B_lo.shape[0]
    x = position_coords(K, read_length) if read_length is not None else np.arange(K)
    n_end = (K - 1) // 2

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Total damage rate
    ax = axes[0]
    r_lo = total_damage_rate(B_lo)
    r_hi = total_damage_rate(B_hi)
    ax.plot(x[:n_end], r_lo[:n_end], 'o-', color='steelblue', label=f"alpha={alpha_lo:.6f} (lower)")
    ax.plot(x[:n_end], r_hi[:n_end], 's--', color='coral', label=f"alpha={alpha_hi:.6f} (upper)")
    ax.plot(x[n_end], r_lo[n_end], 'o', color='steelblue', markersize=8)
    ax.plot(x[n_end], r_hi[n_end], 's', color='coral', markersize=8)
    ax.plot(x[n_end+1:], r_lo[n_end+1:], 'o-', color='steelblue')
    ax.plot(x[n_end+1:], r_hi[n_end+1:], 's--', color='coral')
    ax.set_xlabel("Position in read (bp)")
    ax.set_ylabel("Total damage/error rate")
    ax.set_title(f"{prefix.upper()} — Total damage rate by position")
    ax.legend(fontsize=9)
    ax.set_ylim(bottom=0)

    # C->T damage rate
    ax = axes[1]
    ct_lo = ct_damage_rate(B_lo)
    ct_hi = ct_damage_rate(B_hi)
    ax.plot(x[:n_end], ct_lo[:n_end], 'o-', color='steelblue', label=f"alpha={alpha_lo:.6f} (lower)")
    ax.plot(x[:n_end], ct_hi[:n_end], 's--', color='coral', label=f"alpha={alpha_hi:.6f} (upper)")
    ax.plot(x[n_end], ct_lo[n_end], 'o', color='steelblue', markersize=8)
    ax.plot(x[n_end], ct_hi[n_end], 's', color='coral', markersize=8)
    ax.plot(x[n_end+1:], ct_lo[n_end+1:], 'o-', color='steelblue')
    ax.plot(x[n_end+1:], ct_hi[n_end+1:], 's--', color='coral')
    ax.set_xlabel("Position in read (bp)")
    ax.set_ylabel("C → T damage rate")
    ax.set_title(f"{prefix.upper()} — C→T damage rate by position")
    ax.legend(fontsize=9)
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    path = os.path.join(out_dir, f'{prefix}_damage_profiles.png')
    fig.savefig(path, dpi=150)
    print(f"  Saved {path}")
    plt.close(fig)


def plot_feasibility_scan(Q, C, D, prefix, alpha_lo, alpha_hi, out_dir, alpha_upper=2.0):
    """Plot min B entry as a function of alpha to visualize the feasible interval."""
    C_inv = inv(C)
    alphas = np.linspace(0, alpha_upper, 500)
    min_entries = []
    for a in alphas:
        A = compute_A(a, Q)
        A_inv = inv(A)
        min_val = min(np.min(A_inv @ D[k] @ C_inv) for k in range(D.shape[0]))
        min_entries.append(min_val)
    min_entries = np.array(min_entries)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(alphas, min_entries, 'b-', lw=1.5)
    ax.axhline(0, color='k', ls='-', lw=0.5)
    if alpha_lo is not None:
        ax.axvline(alpha_lo, color='steelblue', ls='--', lw=1, label=f'alpha_lo={alpha_lo:.6f}')
    if alpha_hi is not None:
        ax.axvline(alpha_hi, color='coral', ls='--', lw=1, label=f'alpha_hi={alpha_hi:.6f}')
    ax.set_xlabel('alpha')
    ax.set_ylabel('min entry of B(alpha)')
    ax.set_title(f'{prefix.upper()} — Feasibility scan')
    ax.legend()
    fig.tight_layout()
    path = os.path.join(out_dir, f'{prefix}_feasibility_scan.png')
    fig.savefig(path, dpi=150)
    print(f"  Saved {path}")
    plt.close(fig)


def run_one(prefix, data_dir, out_dir, n_reads=None, read_length=None,
            n_ref=None, alpha_upper=None):
    """Full analysis for one dataset."""
    Q = np.load(os.path.join(data_dir, f'{prefix}_A.npy'))
    C = np.load(os.path.join(data_dir, f'{prefix}_C.npy'))
    D = np.load(os.path.join(data_dir, f'{prefix}_D.npy'))
    K = D.shape[0]

    use_ref_counts = n_ref is not None

    if not use_ref_counts:
        if n_reads is None:
            n_reads = PARAMS.get(prefix, {}).get('n_reads')
        if read_length is None:
            read_length = PARAMS.get(prefix, {}).get('read_length')
        if n_reads is None or read_length is None:
            raise ValueError(
                f"Prefix '{prefix}' not in PARAMS and no --n-reads/--read-length or "
                f"--n-ref provided. Supply one of these on the command line.")
        desc = f"n={n_reads}, L={read_length}"
    else:
        n_ref = np.asarray(n_ref, dtype=float)
        desc = f"n_ref total={n_ref.sum():.0f}"

    print(f"\n{'='*60}")
    print(f"  {prefix.upper()} analysis  (K={K}, {desc})")
    print(f"{'='*60}")

    alpha_max = compute_alpha_max(Q)
    print(f"  alpha_max = 1/max_i(-Q_ii) = {alpha_max:.6f}")

    # --- Step 1: Feasible interval ---
    feas = find_feasible_interval(Q, C, D, alpha_upper=alpha_max, n_grid=10000)
    alpha_lo = feas[0] if feas[0] is not None else 0.0
    alpha_hi = feas[1] if feas[1] is not None else 0.0
    print(f"  Feasible interval: [{alpha_lo:.8f}, {alpha_hi:.8f}]")

    # --- Step 2: Profile likelihood ---
    print(f"\n  Profile likelihood estimation:")
    if use_ref_counts:
        pl_result = estimate_alpha_mll_ref(
            Q, C, D, n_ref=n_ref, alpha_upper=alpha_upper, n_grid=100, verbose=True
        )
    else:
        pl_result = estimate_alpha_mll(
            Q, C, D, n_reads=n_reads, read_length=read_length,
            alpha_upper=alpha_upper, n_grid=100, verbose=True
        )

    plot_profile_likelihood(
        pl_result,
        title=f'{prefix.upper()} — Profile log-likelihood ({desc})',
        save_path=os.path.join(out_dir, f'{prefix}_profile_likelihood.png')
    )

    ci_lo = pl_result.ci_lo
    ci_hi = pl_result.ci_hi

    # --- Step 3: B matrices at 95% CI boundaries ---
    # Inside the feasible interval the analytic inverse A^-1 D C^-1 is non-negative and is the
    # MLE (it reproduces D exactly). Outside it, the analytic B has negative entries, so use
    # the constrained MLL solver instead.
    if use_ref_counts:
        def solve_B(alpha):
            return profile_loglik_ref(alpha, Q, C, D, n_ref)[1]
    else:
        _, n_per_cat = compute_counts(D, n_reads, read_length)
        def solve_B(alpha):
            return profile_loglik(alpha, Q, C, D, n_per_cat)[1]

    B_ci_lo = compute_B_analytic(ci_lo, Q, C, D)
    if B_ci_lo.min() >= -1e-10:
        B_lo_method = 'analytic'
        B_ci_lo = np.maximum(B_ci_lo, 0.0)
        B_ci_lo /= B_ci_lo.sum(axis=2, keepdims=True)
    else:
        B_lo_method = 'constrained'
        print(f"  CI lower bound is outside the feasible interval (analytic B has entries down to "
              f"{B_ci_lo.min():.3g}); solving B_lo with the constrained solver")
        B_ci_lo = solve_B(ci_lo)

    B_ci_hi = solve_B(ci_hi)

    np.save(os.path.join(out_dir, f'{prefix}_B_lo.npy'), B_ci_lo)
    np.save(os.path.join(out_dir, f'{prefix}_B_hi.npy'), B_ci_hi)
    print(f"  Saved {prefix}_B_lo.npy (at CI lower={ci_lo:.6f}, {B_lo_method})")
    print(f"  Saved {prefix}_B_hi.npy (at CI upper={ci_hi:.6f})")

    # Save alpha results
    savez_kwargs = dict(
        alpha_lo=alpha_lo, alpha_hi=alpha_hi,
        alpha_mle=pl_result.alpha_hat,
        ci_lo=ci_lo, ci_hi=ci_hi,
        loglik_max=pl_result.loglik_max,
        profile_alphas=pl_result.profile_alphas,
        profile_logliks=pl_result.profile_logliks,
        alpha_max=alpha_max,
        alpha_upper=pl_result.alpha_upper,
        ci_hi_censored=pl_result.ci_hi_censored,
        B_lo_method=B_lo_method,
    )
    if use_ref_counts:
        savez_kwargs['n_ref'] = n_ref
    else:
        savez_kwargs['n_reads'] = n_reads
        savez_kwargs['read_length'] = read_length
    np.savez(os.path.join(out_dir, f'{prefix}_alpha.npz'), **savez_kwargs)
    print(f"  Saved {prefix}_alpha.npz")

    # --- Step 4: Reconstruction errors ---
    A_ci_lo = compute_A(ci_lo, Q)
    A_ci_hi = compute_A(ci_hi, Q)
    recon_lo = np.array([np.linalg.norm(A_ci_lo @ B_ci_lo[k] @ C - D[k], 'fro') for k in range(K)])
    recon_hi = np.array([np.linalg.norm(A_ci_hi @ B_ci_hi[k] @ C - D[k], 'fro') for k in range(K)])
    print(f"\n  Reconstruction error (CI lower): mean={recon_lo.mean():.2e}, max={recon_lo.max():.2e}")
    print(f"  Reconstruction error (CI upper): mean={recon_hi.mean():.2e}, max={recon_hi.max():.2e}")

    # --- Step 5: Plots ---
    rl_for_plot = read_length if not use_ref_counts else None
    plot_damage_profiles(B_ci_lo, B_ci_hi, prefix, ci_lo, ci_hi, rl_for_plot, out_dir)
    plot_feasibility_scan(Q, C, D, prefix, alpha_lo, alpha_hi, out_dir, alpha_upper=alpha_max)

    # --- Print selected B matrices ---
    nucs = 'ACGT'
    for label, idx, B_mat, alpha_val in [
        ("CI lower", 0, B_ci_lo, ci_lo), ("CI upper", 0, B_ci_hi, ci_hi),
        ("CI lower", K//2, B_ci_lo, ci_lo), ("CI upper", K//2, B_ci_hi, ci_hi)
    ]:
        cat_label = "5' pos 1" if idx == 0 else "middle"
        print(f"\n  B^({cat_label}) at alpha={alpha_val:.8f} ({label}):")
        print(f"       {'  '.join(f'{n:>8s}' for n in nucs)}")
        for i, ni in enumerate(nucs):
            print(f"    {ni}  {'  '.join(f'{B_mat[idx,i,j]:8.5f}' for j in range(4))}")

    return {
        'prefix': prefix,
        'alpha_lo': alpha_lo, 'alpha_hi': alpha_hi,
        'alpha_mle': pl_result.alpha_hat,
        'ci_lo': ci_lo, 'ci_hi': ci_hi,
        'alpha_upper': pl_result.alpha_upper,
        'ci_hi_censored': pl_result.ci_hi_censored,
        'B_lo_method': B_lo_method,
    }


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        description="Run damage estimation analysis.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # existing behaviour (ss and ds with hardcoded params)\n"
            "  python run_analysis.py\n\n"
            "  # new dataset with n_reads / read_length\n"
            "  python run_analysis.py --prefixes mypref --data-dir /path/to/matrices"
            " --n-reads 50000 --read-length 42\n\n"
            "  # new dataset with per-category reference counts (K x 4 .npy)\n"
            "  python run_analysis.py --prefixes mypref --data-dir /path/to/matrices"
            " --n-ref mypref /path/to/mypref_nref.npy\n\n"
            "  # ss and ds each with their own reference count files\n"
            "  python run_analysis.py"
            " --n-ref ss ss_nref.npy --n-ref ds ds_nref.npy\n"
        ),
    )
    parser.add_argument(
        '--prefixes', nargs='+', default=['ss', 'ds'], metavar='PREFIX',
        help="Dataset prefix(es) to analyse (default: ss ds)")
    parser.add_argument(
        '--data-dir', default=script_dir, metavar='DIR',
        help="Directory containing {prefix}_A/C/D.npy files (default: script directory)")
    parser.add_argument(
        '--out-dir', default=os.path.join(script_dir, 'data_analysis_results'), metavar='DIR',
        help="Output directory (default: data_analysis_results/)")
    parser.add_argument(
        '--n-reads', type=int, default=None, metavar='N',
        help="Number of reads; overrides PARAMS default for all prefixes")
    parser.add_argument(
        '--read-length', type=float, default=None, metavar='L',
        help="Read length in bp; overrides PARAMS default for all prefixes")
    parser.add_argument(
        '--n-ref', nargs=2, action='append', default=[], metavar=('PREFIX', 'PATH'),
        help="Reference count file (K x 4 .npy) for a specific prefix. "
             "Can be repeated: --n-ref ss ss_nref.npy --n-ref ds ds_nref.npy")
    parser.add_argument(
        '--alpha-upper', type=float, default=None, metavar='A',
        help="Upper bound for alpha in the profile likelihood scan. Defaults to, and cannot "
             "exceed, alpha_max = 1/max_i(-Q_ii), the largest valid alpha")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    n_ref_map = {}
    for prefix, path in args.n_ref:
        n_ref_map[prefix] = np.load(path)
        print(f"  Loaded n_ref for '{prefix}' from {path}: shape={n_ref_map[prefix].shape}")

    results = {}
    for prefix in args.prefixes:
        results[prefix] = run_one(
            prefix,
            data_dir=args.data_dir,
            out_dir=args.out_dir,
            n_reads=args.n_reads,
            read_length=args.read_length,
            n_ref=n_ref_map.get(prefix),
            alpha_upper=args.alpha_upper,
        )

    print(f"\n{'='*60}")
    print("  Summary")
    print(f"{'='*60}")
    for prefix in args.prefixes:
        r = results[prefix]
        print(f"  {prefix.upper()}:")
        print(f"    Feasible interval: [{r['alpha_lo']:.8f}, {r['alpha_hi']:.8f}]")
        print(f"    MLE: {r['alpha_mle']:.8f}")
        print(f"    95% CI: [{r['ci_lo']:.8f}, {r['ci_hi']:.8f}]")
        print(f"    B_lo solved with: {r['B_lo_method']}")
        if r['ci_hi_censored']:
            print(f"    WARNING: upper CI bound censored at alpha = {r['alpha_upper']:.6g}; "
                  f"B_hi is evaluated there, not at a likelihood-based bound")
    print(f"\nAll results saved to {args.out_dir}")


if __name__ == '__main__':
    main()
