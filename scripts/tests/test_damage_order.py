#!/usr/bin/env python3
"""
Tests that damage estimation (step 8) and the nucleotide likelihoods (step 9) use the
model D^(k) = A(alpha) B^(k) C.

All matrices are row-stochastic, M[from, to], so the product reads left to right in the
physical order reference -> evolution (A) -> damage (B) -> sequencing error (C) -> read.
The fixture uses a large, base-specific sequencing error so that putting C or A on the
wrong side changes the results by much more than the tolerances.

Run with `pytest scripts/tests` or `python scripts/tests/test_damage_order.py`.
"""

import ast
import math
import os
import sys
import tempfile
from types import SimpleNamespace

os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
from scipy.linalg import expm

SCRIPTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SCRIPTS)

from damage_estimator import (compute_B_analytic, find_feasible_interval,
                              solve_inner_wls, solve_inner_mll)
from profile_likelihood import (loglik_multinomial, solve_inner_mll_counts,
                                solve_inner_mll_counts_ref)
import run_analysis

ALPHA = 0.05
N_END = 1                 # bins: 5' position 1, middle, 3' position 1
K = 2 * N_END + 1
PI = np.array([0.32, 0.18, 0.18, 0.32])
A_, C_, G_, T_ = range(4)


def gtr_Q(pi, rates):
    """Normalised GTR rate matrix (as generate_matrices.build_A_matrix); rates AC AG AT CG CT GT."""
    Q = np.zeros((4, 4))
    for (i, j), r in zip([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], rates):
        Q[i, j] = r * pi[j]
        Q[j, i] = r * pi[i]
    np.fill_diagonal(Q, -Q.sum(1))
    return Q / -(pi @ np.diag(Q))


def error_matrix(err):
    """Sequencing error matrix with C[i, i] = 1 - err[i] and C[i, j] = err[i] / 3."""
    C = np.tile(np.asarray(err, dtype=float)[:, None] / 3, (1, 4))
    np.fill_diagonal(C, 1 - np.asarray(err))
    return C


def damage_matrices():
    """Double-stranded damage, with exact zeros for every other substitution."""
    B = np.zeros((K, 4, 4))
    B[0, C_, T_], B[0, G_, T_] = 0.25, 0.03     # 5' end
    B[1, C_, T_], B[1, G_, A_] = 0.01, 0.01     # middle
    B[2, G_, A_], B[2, C_, A_] = 0.25, 0.03     # 3' end
    for k in range(K):
        np.fill_diagonal(B[k], 1 - B[k].sum(1))
    return B


Q = gtr_Q(PI, [1.0, 3.8, 0.8, 1.1, 3.9, 1.0])
C = error_matrix([0.02, 0.015, 0.01, 0.025])
A = expm(ALPHA * Q)
B = damage_matrices()
D = np.array([A @ B[k] @ C for k in range(K)])          # noiseless observations
N_REF = np.round(1e6 * np.tile(PI, (K, 1)))


def test_kron_vectorisation():
    """Row-major vec(A B C) = kron(A, C^T) vec(B), the identity the inner solvers rely on."""
    np.testing.assert_allclose(np.kron(A, C.T) @ B[0].ravel(), (A @ B[0] @ C).ravel(), atol=1e-15)


def test_monte_carlo_physical_order():
    """Base-by-base draws through evolution, damage and sequencing error give A B C, not C B A."""
    rng = np.random.default_rng(1)
    n = 200_000
    alpha = 0.3                                # large, so the two orders are far apart
    A_mc = expm(alpha * Q)

    def draw(M, states):
        return (M[states].cumsum(1) < rng.random(len(states))[:, None]).sum(1).clip(0, 3)

    for i in range(4):
        x = draw(C, draw(B[0], draw(A_mc, np.full(n, i))))
        freq = np.bincount(x, minlength=4) / n
        for P, should_match in ((A_mc @ B[0] @ C, True), (C @ B[0] @ A_mc, False)):
            z = np.abs(freq - P[i]) / np.sqrt(P[i] * (1 - P[i]) / n)
            assert (z.max() < 5) == should_match, (i, should_match, z)


def test_analytic_B_recovers_truth():
    np.testing.assert_allclose(compute_B_analytic(ALPHA, Q, C, D), B, atol=1e-12)


def test_inner_solvers_recover_truth():
    for k in range(K):
        for solve in (solve_inner_wls, solve_inner_mll):
            np.testing.assert_allclose(solve(D[k], C, A)[0], B[k], atol=1e-6)
        np.testing.assert_allclose(
            solve_inner_mll_counts(D[k], C, A, N_REF[k].sum(), base_freq=PI)[0], B[k], atol=1e-6)
        np.testing.assert_allclose(solve_inner_mll_counts_ref(D[k], C, A, N_REF[k])[0], B[k], atol=1e-6)


def test_loglik_multinomial_model():
    """loglik_multinomial uses P = A B C: at the true B it equals the maximum, sum n f D log D."""
    n_per_cat = N_REF.sum(1)
    expected = sum(n_per_cat[k] * PI[i] * np.sum(D[k, i] * np.log(D[k, i]))
                   for k in range(K) for i in range(4))
    np.testing.assert_allclose(loglik_multinomial(D, C, B, A, n_per_cat, base_freq=PI),
                               expected, rtol=1e-12)


def test_feasible_interval():
    """With exact zeros in B, the analytic B leaves the feasible set exactly at the true alpha."""
    lo, hi = find_feasible_interval(Q, C, D, alpha_upper=1.0, n_grid=2000)
    assert lo == 0.0
    np.testing.assert_allclose(hi, ALPHA, rtol=1e-6)


def test_run_analysis_end_to_end():
    """run_analysis.run_one (--n-ref path) brackets the true alpha and recovers B at ci_hi."""
    with tempfile.TemporaryDirectory() as tmp:
        for name, M in (('A', Q), ('C', C), ('D', D)):
            np.save(os.path.join(tmp, f'sim_{name}.npy'), M)
        res = run_analysis.run_one('sim', tmp, tmp, n_ref=N_REF)
        B_lo = np.load(os.path.join(tmp, 'sim_B_lo.npy'))
        B_hi = np.load(os.path.join(tmp, 'sim_B_hi.npy'))

    np.testing.assert_allclose(res['alpha_hi'], ALPHA, rtol=1e-6)     # feasible interval
    assert res['ci_lo'] <= ALPHA <= res['ci_hi'] <= 1.1 * ALPHA
    assert res['B_lo_method'] == 'analytic'
    # At ci_lo = 0, A = I and B_lo = D C^-1 (= A B): only the position of C is tested here
    assert res['ci_lo'] == 0.0
    np.testing.assert_allclose(B_lo, D @ np.linalg.inv(C), atol=1e-10)
    np.testing.assert_allclose(B_hi, B, atol=5e-3)


def _load_step9_functions():
    """Load the step-9 likelihood code without running the script (it parses arguments at import)."""
    path = os.path.join(SCRIPTS, 'generateRatePlacerInputs_msaReads_damageMatrix.py')
    with open(path) as f:
        tree = ast.parse(f.read())
    wanted = {'generate_error_profile', 'quality_to_error_prob', 'build_positional_matrices_from_npy'}
    funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
    assert {n.name for n in funcs} == wanted
    ns = {'math': math, 'np': np,
          'args': SimpleNamespace(no_error_profile=False, quality_score_only=False,
                                  positional_matrix_bases=N_END)}
    exec(compile(ast.Module(body=funcs, type_ignores=[]), path, 'exec'), ns)
    return ns


def test_step9_likelihoods_match_estimation_model():
    """Step 9's P(X | Y) (Eq. 7), averaged over Y ~ A(alpha)[ref], is the step-8 model (A B C)[ref, X]."""
    ns = _load_step9_functions()
    q = 17                                    # Phred 17: e = 0.020, the same for every base
    e = 10 ** (-q / 10)
    C_eq = error_matrix([e] * 4)
    D_eq = np.array([A @ B[k] @ C_eq for k in range(K)])
    positional = ns['build_positional_matrices_from_npy'](B, N_END)

    bins = [0, 1, 1, 1, 2]                    # read positions 0 / 1-3 / 4 -> 5' / middle / 3'
    for obs in 'ACGT':
        seq = [obs] * len(bins)
        entries = ns['generate_error_profile'](
            seq, chr(q + 33) * len(seq), list(range(len(seq))), 99, 'phred33',
            None, B[N_END], 'A' * len(seq), None, positional, False, positional_bases=N_END)
        assert len(entries) == len(seq)
        x = 'ACGT'.index(obs)
        for p, line in enumerate(entries):
            L = np.exp([float(v) for v in line.split('\t')[2:]])     # P(X = obs | Y), Y = A C G T
            k = bins[p]
            # Eq. 7 directly, then the estimator's own likelihood of reference i -> observed X
            np.testing.assert_allclose(L, B[k] @ C_eq[:, x], rtol=1e-5)
            for i in range(4):
                onehot = np.zeros((K, 4, 4))
                onehot[k, i, x] = 1.0
                weight = np.zeros(K)
                weight[k] = 1.0
                base_freq = np.eye(4)[i]
                P_model = np.exp(loglik_multinomial(onehot, C_eq, B, A, weight, base_freq=base_freq))
                np.testing.assert_allclose(A[i] @ L, P_model, rtol=1e-5)
                np.testing.assert_allclose(P_model, D_eq[k, i, x], rtol=1e-12)


if __name__ == '__main__':
    tests = [(name, f) for name, f in sorted(globals().items()) if name.startswith('test_')]
    failed = 0
    for name, f in tests:
        try:
            f()
            print(f'PASS  {name}')
        except Exception as exc:          # report every test, not just the first failure
            failed += 1
            print(f'FAIL  {name}: {type(exc).__name__}: {exc}')
    print(f'\n{len(tests) - failed}/{len(tests)} passed')
    sys.exit(1 if failed else 0)
