"""Numerical-correctness tests for the spinner iteration."""

from __future__ import annotations

import numpy as np
import pytest

from winner.core import (
    initial_score_from_adj,
    spinner_iteration,
    spinner_iteration_batch,
)


def _toy_net():
    # 4-node weighted, undirected graph.
    A = np.array(
        [
            [0.0, 0.7, 0.0, 0.3],
            [0.7, 0.0, 0.5, 0.0],
            [0.0, 0.5, 0.0, 0.8],
            [0.3, 0.0, 0.8, 0.0],
        ]
    )
    return A


def test_initial_score_matches_wdeg2_over_deg():
    A = _toy_net()
    v0 = initial_score_from_adj(A)
    wdeg = A.sum(axis=0)
    deg = np.sign(A).sum(axis=0)
    expected = wdeg**2 / deg
    np.testing.assert_allclose(v0, expected, rtol=1e-12)


def test_initial_score_handles_isolated_nodes():
    A = np.zeros((3, 3))
    A[0, 1] = A[1, 0] = 0.5
    v0 = initial_score_from_adj(A)
    assert v0[2] == 0.0
    assert v0[0] > 0.0


def test_spinner_converges_to_fixed_point():
    A = _toy_net()
    v0 = initial_score_from_adj(A)
    v = spinner_iteration(A, v0, max_iter=500)

    row_sum = A.sum(axis=1, keepdims=True)
    A_norm = A / row_sum
    sigma = 0.85
    rhs = (1 - sigma) * v0 + sigma * A_norm.T @ v
    np.testing.assert_allclose(v, rhs, rtol=1e-10, atol=1e-12)


def test_batch_matches_single():
    A = _toy_net()
    v0 = initial_score_from_adj(A)
    single = spinner_iteration(A, v0)
    batched = spinner_iteration_batch(
        np.stack([A, A, A], axis=0),
        np.stack([v0, v0, v0], axis=0),
    )
    for b in range(3):
        np.testing.assert_allclose(batched[b], single, rtol=1e-12)


def test_torch_batch_matches_numpy_when_available():
    torch = pytest.importorskip("torch")
    A = _toy_net()
    v0 = initial_score_from_adj(A)
    single = spinner_iteration(A, v0)
    from winner.core import spinner_iteration_torch_batch

    batched = spinner_iteration_torch_batch(
        np.stack([A, A], axis=0),
        np.stack([v0, v0], axis=0),
        device="cpu",
        dtype="float64",
    )
    np.testing.assert_allclose(batched[0], single, rtol=1e-10)


def test_sparse_batch_matches_dense():
    from winner.core import spinner_iteration_sparse_batch

    A = _toy_net()
    v0 = initial_score_from_adj(A)
    dense = spinner_iteration_batch(
        np.stack([A, A, A], axis=0), np.stack([v0, v0, v0], axis=0)
    )
    sparse = spinner_iteration_sparse_batch(
        np.stack([A, A, A], axis=0), np.stack([v0, v0, v0], axis=0)
    )
    np.testing.assert_allclose(sparse, dense, rtol=1e-10)


def test_dispatcher_autoselects_sparse_for_sparse_input():
    from winner.core import spinner_batch

    # Tiny, sparse network — density << 5%, CPU device → should pick sparse path.
    N = 50
    A = np.zeros((N, N))
    rng = np.random.default_rng(0)
    ii = rng.integers(0, N, size=10)
    jj = rng.integers(0, N, size=10)
    for i, j in zip(ii, jj):
        if i != j:
            w = rng.uniform(0.1, 1.0)
            A[i, j] = w
            A[j, i] = w
    v0 = initial_score_from_adj(A)
    B = 5
    out = spinner_batch(
        np.stack([A] * B, axis=0),
        np.stack([v0] * B, axis=0),
        device="cpu",
    )
    assert out.shape == (B, N)
    assert np.all(np.isfinite(out))
