"""Sanity tests for the degree-preserving rewire."""

from __future__ import annotations

import numpy as np

from winner.randomize import generate_random_network, generate_random_networks_batch


def _ring_net(n=10):
    A = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        A[i, (i + 1) % n] = 1.0
        A[(i + 1) % n, i] = 1.0
    return A


def test_preserves_degree_sequence():
    A = _ring_net(20)
    rng_seed = 123
    R = generate_random_network(A, seed=rng_seed)
    np.testing.assert_array_equal(
        np.sign(A).sum(axis=0), np.sign(R).sum(axis=0)
    )


def test_preserves_edge_count():
    A = _ring_net(30)
    R = generate_random_network(A, seed=42)
    n_edges_A = int((np.triu(A, 1) > 0).sum())
    n_edges_R = int((np.triu(R, 1) > 0).sum())
    assert n_edges_A == n_edges_R


def test_weighted_rewire_samples_from_pool():
    A = _ring_net(12)
    rng = np.random.default_rng(0)
    weights = rng.uniform(0.1, 1.0, size=int((np.triu(A, 1) > 0).sum()))
    triu_idx = np.where(np.triu(A, 1) > 0)
    weighted = A.copy()
    weighted[triu_idx] = weights
    weighted[triu_idx[1], triu_idx[0]] = weights
    R = generate_random_network(weighted, edge_weights=weights, seed=7)
    observed = R[np.triu(R, 1) > 0]
    # each observed weight is drawn from the pool
    assert set(np.round(observed, 10).tolist()).issubset(
        set(np.round(weights, 10).tolist())
    )


def test_batch_shape_and_determinism():
    A = _ring_net(8)
    stack_a = generate_random_networks_batch(A, n=5, seed=1, n_jobs=1)
    stack_b = generate_random_networks_batch(A, n=5, seed=1, n_jobs=1)
    assert stack_a.shape == (5, 8, 8)
    np.testing.assert_array_equal(stack_a, stack_b)
