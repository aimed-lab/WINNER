"""Degree-preserving random-network generation.

Ports ``sym_generate_srand.m``. Given a binary symmetric adjacency, repeatedly
swap pairs of edges (``(v1,v2)`` and ``(v3,v4)`` become ``(v1,v3)`` and
``(v2,v4)``) while ensuring endpoints stay distinct and no parallel edge is
created. The resulting graph has the same degree sequence as the input.

If :mod:`numba` is installed we JIT the hot loop for a ~30x speed-up; otherwise
we fall back to pure NumPy (still usable, just slower).
"""

from __future__ import annotations

from typing import Optional

import numpy as np

try:
    from numba import njit
    _HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    _HAVE_NUMBA = False

    def njit(*args, **kwargs):
        def wrap(f):
            return f
        if args and callable(args[0]):
            return args[0]
        return wrap


@njit(cache=True)
def _rewire_binary(
    adj: np.ndarray,
    edges_i: np.ndarray,
    edges_j: np.ndarray,
    ntry: int,
    seed: int,
) -> np.ndarray:
    """In-place double-edge swap on a binary symmetric matrix.

    ``edges_i[k], edges_j[k]`` track the upper-triangle endpoint of each edge.
    The array is mutated to reflect accepted swaps. Returns ``adj``.
    """
    np.random.seed(seed)
    ne = edges_i.shape[0]
    for _ in range(ntry):
        e1 = int(np.random.randint(0, ne))
        e2 = int(np.random.randint(0, ne))
        v1 = edges_i[e1]
        v2 = edges_j[e1]
        v3 = edges_i[e2]
        v4 = edges_j[e2]
        if v1 == v3 or v1 == v4 or v2 == v3 or v2 == v4:
            continue
        if np.random.random() > 0.5:
            # try edges (v1,v3) and (v2,v4)
            if adj[v1, v3] == 0 and adj[v2, v4] == 0:
                adj[v1, v2] = 0
                adj[v2, v1] = 0
                adj[v3, v4] = 0
                adj[v4, v3] = 0
                adj[v1, v3] = 1
                adj[v3, v1] = 1
                adj[v2, v4] = 1
                adj[v4, v2] = 1
                edges_i[e1] = v1
                edges_j[e1] = v3
                edges_i[e2] = v2
                edges_j[e2] = v4
        else:
            # try edges (v1,v4) and (v2,v3)
            if adj[v1, v4] == 0 and adj[v2, v3] == 0:
                adj[v1, v2] = 0
                adj[v2, v1] = 0
                adj[v3, v4] = 0
                adj[v4, v3] = 0
                adj[v1, v4] = 1
                adj[v4, v1] = 1
                adj[v2, v3] = 1
                adj[v3, v2] = 1
                edges_i[e1] = v1
                edges_j[e1] = v4
                edges_i[e2] = v2
                edges_j[e2] = v3
    return adj


def generate_random_network(
    adj: np.ndarray,
    ntry: Optional[int] = None,
    edge_weights: Optional[np.ndarray] = None,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Return one degree-preserved random network.

    If ``edge_weights`` is given, each edge of the rewired network is assigned
    a weight sampled uniformly with replacement from that pool (matches the
    MATLAB p-value routine).
    """
    rng = np.random.default_rng(seed)
    binary = (np.asarray(adj) != 0).astype(np.int8)
    ii, jj = np.where(np.triu(binary, k=1) > 0)
    edges_i = ii.astype(np.int64)
    edges_j = jj.astype(np.int64)
    if edges_i.size == 0:
        return np.zeros_like(adj, dtype=np.float64)

    n_edges = edges_i.size
    tries = int(ntry) if ntry is not None else 4 * n_edges
    jit_seed = int(rng.integers(0, 2**31 - 1))
    rewired = _rewire_binary(binary.copy(), edges_i.copy(), edges_j.copy(), tries, jit_seed)

    if edge_weights is None:
        return rewired.astype(np.float64)

    pool = np.asarray(edge_weights, dtype=np.float64)
    weights = rng.choice(pool, size=n_edges, replace=True)
    ii2, jj2 = np.where(np.triu(rewired, k=1) > 0)
    out = np.zeros_like(rewired, dtype=np.float64)
    out[ii2, jj2] = weights
    out[jj2, ii2] = weights
    return out


def generate_random_networks_batch(
    adj: np.ndarray,
    n: int,
    edge_weights: Optional[np.ndarray] = None,
    ntry: Optional[int] = None,
    seed: Optional[int] = None,
    n_jobs: int = 1,
    backend: str = "threading",
) -> np.ndarray:
    """Generate ``n`` random networks, optionally in parallel.

    Returns an array of shape ``(n, N, N)``.

    Threading backend is default: the inner rewire loop is Numba-JIT which
    releases the GIL, and threads avoid the pickling / worker-spawn cost
    that dominates process pools on small-to-medium networks.
    """
    rng = np.random.default_rng(seed)
    child_seeds = rng.integers(0, 2**31 - 1, size=n).tolist()

    def _one(s: int) -> np.ndarray:
        return generate_random_network(adj, ntry=ntry, edge_weights=edge_weights, seed=int(s))

    if n_jobs == 1 or n < 4:
        stack = np.stack([_one(s) for s in child_seeds], axis=0)
        return stack

    from joblib import Parallel, delayed

    results = Parallel(n_jobs=n_jobs, backend=backend)(
        delayed(_one)(s) for s in child_seeds
    )
    return np.stack(results, axis=0)
