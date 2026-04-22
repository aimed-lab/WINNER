"""Core spinner / personalized-PageRank iteration.

Ports ``spinnerIteration.m`` from the MATLAB WINNER implementation. Given a
weighted, undirected adjacency matrix ``adj`` and an initial score vector
``v0``, the spinner score is the fixed-point of

    v_{t+1} = (1 - sigma) * v_0 + sigma * A^T v_t

where ``A`` is row-normalised (rows that sum to zero are left untouched).

We expose five entry points, auto-selected by :func:`spinner_batch` based on
network density and device availability:

* :func:`spinner_iteration` — NumPy, one network.
* :func:`spinner_iteration_batch` — NumPy dense, ``(B, N, N)`` via ``matmul``
  (BLAS-backed, unlike ``einsum``).
* :func:`spinner_iteration_sparse_batch` — SciPy CSR per network (winning path
  for PPI networks with < ~5% density).
* :func:`spinner_iteration_torch_batch` — PyTorch dense ``bmm`` on CPU /
  CUDA / MPS (float32 by default off-CPU).
* :func:`spinner_iteration_torch_sparse_batch` — PyTorch sparse CSR / COO
  BMM on GPU when the network is sparse AND the device supports it.

All five are numerically equivalent up to floating-point order of operations.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .backend import Device, resolve_device


DEFAULT_MAX_ITER = 100
DEFAULT_SIGMA = 0.85
SPARSE_DENSITY_THRESHOLD = 0.05  # < 5% density → sparse path wins


def _row_normalize(adj: np.ndarray) -> np.ndarray:
    """Divide each row by its sum, leaving zero-sum rows as zero."""
    row_sum = adj.sum(axis=1, keepdims=True)
    safe = np.where(row_sum > 0, row_sum, 1.0)
    out = adj / safe
    out[np.squeeze(row_sum <= 0, axis=1)] = 0.0
    return out


def _row_normalize_batch(adj: np.ndarray) -> np.ndarray:
    row_sum = adj.sum(axis=2, keepdims=True)
    safe = np.where(row_sum > 0, row_sum, 1.0)
    out = adj / safe
    mask = np.squeeze(row_sum <= 0, axis=2)
    out[mask] = 0.0
    return out


def initial_score_from_adj(adj: np.ndarray) -> np.ndarray:
    """Compute ``wdeg^2 / deg`` initial score as in ``RunWinner.m``.

    Matches ``initialScore = exp(2*log(nodeWDeg) - log(nodeDeg))`` with NaN→0.
    """
    w_deg = adj.sum(axis=0)
    deg = np.sign(adj).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        score = np.where(deg > 0, (w_deg ** 2) / deg, 0.0)
    score = np.nan_to_num(score, nan=0.0, posinf=0.0, neginf=0.0)
    return score.astype(np.float64, copy=False)


def initial_score_from_adj_batch(adj_batch: np.ndarray) -> np.ndarray:
    """Batched ``wdeg^2 / deg``. Input ``(B, N, N)`` → output ``(B, N)``.

    We intentionally loop over ``B``: for typical WINNER inputs (``N ≲ 1000``)
    each ``(N, N)`` slice fits in L2/L3 cache and NumPy per-call overhead is
    negligible, whereas vectorising over ``B`` forces a large bool/sign
    temporary that thrashes the page cache.
    """
    B, N, _ = adj_batch.shape
    out = np.empty((B, N), dtype=np.float64)
    for b in range(B):
        out[b] = initial_score_from_adj(adj_batch[b])
    return out


def spinner_iteration(
    adj: np.ndarray,
    initial_score: np.ndarray,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
) -> np.ndarray:
    """NumPy spinner for one network. Returns the final score vector."""
    if adj.ndim != 2 or adj.shape[0] != adj.shape[1]:
        raise ValueError("adj must be square 2-D")
    if initial_score.shape != (adj.shape[0],):
        raise ValueError("initial_score length must equal adj size")

    A = _row_normalize(np.asarray(adj, dtype=np.float64))
    AT = A.T.copy()
    v0 = np.asarray(initial_score, dtype=np.float64)
    v = v0.copy()
    one_minus_sigma_v0 = (1.0 - sigma) * v0
    for _ in range(max_iter - 1):
        v = one_minus_sigma_v0 + sigma * (AT @ v)
    return v


def spinner_iteration_batch(
    adj_batch: np.ndarray,
    initial_batch: np.ndarray,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
) -> np.ndarray:
    """NumPy batched dense spinner. ``adj_batch`` is ``(B, N, N)``.

    Uses ``np.matmul`` (BLAS-backed) rather than ``einsum`` so NumPy dispatches
    to batched gemv via the system's linear-algebra library.
    """
    if adj_batch.ndim != 3:
        raise ValueError("adj_batch must be 3-D (B, N, N)")
    B, N, _ = adj_batch.shape
    if initial_batch.shape != (B, N):
        raise ValueError("initial_batch must be (B, N)")

    A = _row_normalize_batch(np.asarray(adj_batch, dtype=np.float64))
    AT = np.transpose(A, (0, 2, 1))
    v0 = np.asarray(initial_batch, dtype=np.float64)
    v = v0[:, :, None].copy()  # (B, N, 1)
    scaled_v0 = (1.0 - sigma) * v0[:, :, None]
    for _ in range(max_iter - 1):
        v = scaled_v0 + sigma * np.matmul(AT, v)
    return v[:, :, 0]


def spinner_iteration_sparse_batch(
    adj_batch: np.ndarray,
    initial_batch: np.ndarray,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
    n_jobs: int = 1,
) -> np.ndarray:
    """SciPy CSR batched spinner — the fast path for sparse PPI networks.

    Each network is compiled into a CSR matrix once, then multiplied by the
    score vector 99 times. For a 0.2%-dense network this is ~500× fewer FLOPs
    than the dense batched matmul. ``n_jobs > 1`` runs the outer B loop
    across threads (the SciPy kernel releases the GIL).
    """
    from scipy.sparse import csr_matrix

    if adj_batch.ndim != 3:
        raise ValueError("adj_batch must be 3-D (B, N, N)")
    B, N, _ = adj_batch.shape
    out = np.empty((B, N), dtype=np.float64)

    def _one(b: int) -> None:
        A = _row_normalize(np.asarray(adj_batch[b], dtype=np.float64))
        AT = csr_matrix(A.T)
        v0 = np.asarray(initial_batch[b], dtype=np.float64)
        v = v0.copy()
        scaled = (1.0 - sigma) * v0
        for _ in range(max_iter - 1):
            v = scaled + sigma * (AT @ v)
        out[b] = v

    if n_jobs == 1 or B < 4:
        for b in range(B):
            _one(b)
    else:
        from joblib import Parallel, delayed

        Parallel(n_jobs=n_jobs, backend="threading")(
            delayed(_one)(b) for b in range(B)
        )
    return out


def spinner_iteration_torch_batch(
    adj_batch: np.ndarray,
    initial_batch: np.ndarray,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
    device: Device = "auto",
    dtype: Optional[str] = None,
):
    """PyTorch batched dense spinner. Works on CPU, CUDA, and Apple MPS.

    Default precision is float32 everywhere (including CPU) — MKL/oneDNN
    batched ``bmm`` is much faster in float32 and the 100-iter personalized
    PageRank is numerically benign enough for it. Pass ``dtype="float64"``
    for exact MATLAB parity.
    """
    import torch

    resolved = resolve_device(device)
    if dtype is None:
        dtype = "float32"
    torch_dtype = getattr(torch, dtype)
    dev = torch.device(resolved)

    adj_t = torch.as_tensor(adj_batch, dtype=torch_dtype, device=dev)
    v0 = torch.as_tensor(initial_batch, dtype=torch_dtype, device=dev)

    row_sum = adj_t.sum(dim=2, keepdim=True)
    safe = torch.where(row_sum > 0, row_sum, torch.ones_like(row_sum))
    A = adj_t / safe
    A = torch.where(row_sum.expand_as(A) > 0, A, torch.zeros_like(A))
    AT = A.transpose(1, 2).contiguous()

    v = v0.unsqueeze(-1).clone()
    scaled_v0 = (1.0 - sigma) * v0.unsqueeze(-1)
    for _ in range(max_iter - 1):
        v = scaled_v0 + sigma * torch.bmm(AT, v)

    return v.squeeze(-1).detach().cpu().numpy().astype(np.float64, copy=False)


def spinner_iteration_torch_sparse_batch(
    adj_batch: np.ndarray,
    initial_batch: np.ndarray,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
    device: Device = "auto",
    dtype: Optional[str] = None,
):
    """PyTorch sparse batched spinner.

    Builds one big block-diagonal sparse CSR (``B·N × B·N``) and multiplies
    it by a stacked vector ``(B·N,)``. This avoids per-network kernel
    launches on GPU — empirically the winning GPU path for sparse PPI graphs.
    Falls through to per-network sparse matmul on MPS (which does not yet
    support ``torch.sparse.mm`` on block-diag inputs).
    """
    import torch

    resolved = resolve_device(device)
    if dtype is None:
        dtype = "float32"
    torch_dtype = getattr(torch, dtype)
    dev = torch.device(resolved)

    B, N, _ = adj_batch.shape

    # Row-normalise per network (dense, one-shot).
    a = np.asarray(adj_batch, dtype=np.float64)
    row_sum = a.sum(axis=2, keepdims=True)
    safe = np.where(row_sum > 0, row_sum, 1.0)
    A_norm = a / safe
    A_norm[np.squeeze(row_sum <= 0, axis=2)] = 0.0

    # Collect non-zero entries of each AT and offset row/col indices so that
    # block b sits at rows/cols [b*N, (b+1)*N).
    nz_idx = np.nonzero(A_norm)
    b_idx, i_idx, j_idx = nz_idx  # A_norm[b, i, j] -> AT has row=j, col=i
    values = A_norm[b_idx, i_idx, j_idx]
    rows = (b_idx * N + j_idx).astype(np.int64)
    cols = (b_idx * N + i_idx).astype(np.int64)

    indices = torch.from_numpy(np.stack([rows, cols], axis=0)).to(dev)
    vals = torch.as_tensor(values, dtype=torch_dtype, device=dev)
    AT_big = torch.sparse_coo_tensor(indices, vals, (B * N, B * N)).coalesce()

    v0 = torch.as_tensor(initial_batch, dtype=torch_dtype, device=dev).reshape(B * N)
    v = v0.clone()
    scaled_v0 = (1.0 - sigma) * v0
    for _ in range(max_iter - 1):
        v = scaled_v0 + sigma * torch.sparse.mm(AT_big, v.unsqueeze(-1)).squeeze(-1)

    return v.reshape(B, N).detach().cpu().numpy().astype(np.float64, copy=False)


def _density(adj_batch: np.ndarray) -> float:
    B, N, _ = adj_batch.shape
    nz = np.count_nonzero(adj_batch)
    return nz / (B * N * N) if B * N * N else 0.0


def spinner_batch(
    adj_batch: np.ndarray,
    initial_batch: np.ndarray,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
    device: Device = "auto",
    n_jobs: int = 1,
    sparse_threshold: float = SPARSE_DENSITY_THRESHOLD,
    force_dense: bool = False,
    force_sparse: bool = False,
) -> np.ndarray:
    """Auto-select CPU/GPU and sparse/dense path for the batched spinner.

    Rules (unless overridden by ``force_sparse`` / ``force_dense``):

    * ``device == 'cpu'`` and density < ``sparse_threshold`` → SciPy CSR path.
    * ``device in ('cuda', 'mps')`` and density < ``sparse_threshold`` →
      PyTorch sparse block-diag path.
    * Otherwise → dense path (NumPy matmul on CPU, PyTorch bmm on GPU).
    """
    resolved = resolve_device(device)
    is_sparse = force_sparse or (
        not force_dense and _density(adj_batch) < sparse_threshold
    )

    if resolved == "cpu":
        if is_sparse:
            return spinner_iteration_sparse_batch(
                adj_batch, initial_batch, max_iter, sigma, n_jobs=n_jobs
            )
        try:
            return spinner_iteration_torch_batch(
                adj_batch, initial_batch, max_iter, sigma, device="cpu"
            )
        except ImportError:
            return spinner_iteration_batch(
                adj_batch, initial_batch, max_iter, sigma
            )

    # GPU device
    if is_sparse:
        try:
            return spinner_iteration_torch_sparse_batch(
                adj_batch, initial_batch, max_iter, sigma, device=resolved
            )
        except (ImportError, RuntimeError, NotImplementedError):
            pass  # fall back to dense GPU
    return spinner_iteration_torch_batch(
        adj_batch, initial_batch, max_iter, sigma, device=resolved
    )
