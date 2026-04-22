"""Core spinner / personalized-PageRank iteration.

Ports ``spinnerIteration.m`` from the MATLAB WINNER implementation. Given a
weighted, undirected adjacency matrix ``adj`` and an initial score vector
``v0``, the spinner score is the fixed-point of

    v_{t+1} = (1 - sigma) * v_0 + sigma * A^T v_t

where ``A`` is row-normalised (rows that sum to zero are left untouched).

We expose three entry points:

* :func:`spinner_iteration` — NumPy, one network.
* :func:`spinner_iteration_batch` — NumPy, B networks stacked as ``(B, N, N)``.
* :func:`spinner_iteration_torch_batch` — PyTorch, same shape, on any device.

All three are numerically equivalent up to floating-point order of operations.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .backend import Device, resolve_device


DEFAULT_MAX_ITER = 100
DEFAULT_SIGMA = 0.85


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
    """NumPy batched spinner. ``adj_batch`` is ``(B, N, N)``.

    Returns scores of shape ``(B, N)``. Used by the CPU-parallel p-value path
    when PyTorch is not available, or for small batches where GPU overhead
    dominates.
    """
    if adj_batch.ndim != 3:
        raise ValueError("adj_batch must be 3-D (B, N, N)")
    B, N, _ = adj_batch.shape
    if initial_batch.shape != (B, N):
        raise ValueError("initial_batch must be (B, N)")

    A = _row_normalize_batch(np.asarray(adj_batch, dtype=np.float64))
    AT = np.transpose(A, (0, 2, 1))
    v0 = np.asarray(initial_batch, dtype=np.float64)
    v = v0.copy()
    one_minus_sigma_v0 = (1.0 - sigma) * v0
    for _ in range(max_iter - 1):
        v = one_minus_sigma_v0 + sigma * np.einsum("bij,bj->bi", AT, v)
    return v


def spinner_iteration_torch_batch(
    adj_batch: np.ndarray,
    initial_batch: np.ndarray,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
    device: Device = "auto",
    dtype: Optional[str] = None,
):
    """PyTorch batched spinner, device is chosen via :func:`resolve_device`.

    ``dtype`` may be ``"float32"`` or ``"float64"``. Defaults to float32 on
    CUDA/MPS (much faster), float64 on CPU.
    """
    import torch

    resolved = resolve_device(device)
    if dtype is None:
        dtype = "float64" if resolved == "cpu" else "float32"
    torch_dtype = getattr(torch, dtype)
    dev = torch.device(resolved)

    adj_t = torch.as_tensor(adj_batch, dtype=torch_dtype, device=dev)
    v0 = torch.as_tensor(initial_batch, dtype=torch_dtype, device=dev)

    row_sum = adj_t.sum(dim=2, keepdim=True)
    safe = torch.where(row_sum > 0, row_sum, torch.ones_like(row_sum))
    A = adj_t / safe
    A = torch.where(
        row_sum.expand_as(A) > 0, A, torch.zeros_like(A)
    )
    AT = A.transpose(1, 2).contiguous()

    v = v0.clone()
    scaled_v0 = (1.0 - sigma) * v0
    for _ in range(max_iter - 1):
        v = scaled_v0 + sigma * torch.bmm(AT, v.unsqueeze(-1)).squeeze(-1)

    return v.detach().cpu().numpy().astype(np.float64, copy=False)
