"""High-level pipelines: ``run_winner`` and ``run_winner_with_pvalue``.

These mirror ``RunWinner.m`` and ``RunWinner_withPValue.m`` from the reference
MATLAB implementation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .backend import Device, available_devices, resolve_device
from .core import (
    DEFAULT_MAX_ITER,
    DEFAULT_SIGMA,
    initial_score_from_adj,
    initial_score_from_adj_batch,
    spinner_batch,
    spinner_iteration,
)
from .io import build_adjacency
from .pvalue import bh_fdr, expansion_pvalue, ranking_pvalue
from .randomize import generate_random_networks_batch


@dataclass
class WinnerResult:
    gene_names: list[str]
    seed_or_expand: list[str]
    winner_score: np.ndarray

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "geneName": self.gene_names,
                "seedOrExpand": self.seed_or_expand,
                "winnerScore": self.winner_score,
            }
        )


@dataclass
class WinnerPValueResult:
    gene_names: list[str]
    seed_or_expand: list[str]
    winner_score: np.ndarray
    expansion_pval: np.ndarray
    ranking_pval: np.ndarray
    iterations: list[dict] = field(default_factory=list)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "finalGeneList": self.gene_names,
                "finalScore": self.winner_score,
                "expansionPVal": self.expansion_pval,
                "rankingPVal": self.ranking_pval,
            }
        )


def run_winner(
    genes: pd.DataFrame,
    interactions: pd.DataFrame,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
) -> WinnerResult:
    """Simple WINNER (no p-values). Ports ``RunWinner.m``."""
    gene_names = genes["gene"].tolist()
    seed_or_expand = genes["seed_or_expand"].tolist()
    adj = build_adjacency(gene_names, interactions)
    v0 = initial_score_from_adj(adj)
    score = spinner_iteration(adj, v0, max_iter=max_iter, sigma=sigma)
    return WinnerResult(
        gene_names=gene_names,
        seed_or_expand=seed_or_expand,
        winner_score=score,
    )


def _batched_spinner(
    adj_stack: np.ndarray,
    v0_stack: np.ndarray,
    max_iter: int,
    sigma: float,
    device: str,
    chunk: int,
    n_jobs: int = 1,
) -> np.ndarray:
    B = adj_stack.shape[0]
    out = np.empty((B, adj_stack.shape[1]), dtype=np.float64)
    for start in range(0, B, chunk):
        end = min(B, start + chunk)
        out[start:end] = spinner_batch(
            adj_stack[start:end],
            v0_stack[start:end],
            max_iter=max_iter,
            sigma=sigma,
            device=device,
            n_jobs=n_jobs,
        )
    return out


def run_winner_with_pvalue(
    genes: pd.DataFrame,
    interactions: pd.DataFrame,
    global_degree: dict[str, int],
    *,
    total_connected_genes: int = 9967,
    max_expansions: int = 50,
    expansion_alpha: float = 0.05,
    num_random: int = 10000,
    random_seed: Optional[int] = None,
    max_iter: int = DEFAULT_MAX_ITER,
    sigma: float = DEFAULT_SIGMA,
    device: Device = "auto",
    n_jobs: int = -1,
    chunk: int = 500,
    progress: bool = False,
) -> WinnerPValueResult:
    """Full WINNER pipeline with p-values. Ports ``RunWinner_withPValue.m``.

    ``num_random`` controls the empirical null size (MATLAB uses 10000).
    ``n_jobs`` is passed to joblib (``-1`` = all cores) for random-network
    generation. The batched spinner over the null networks uses PyTorch when
    available on ``cuda`` / ``mps``, otherwise a NumPy batched fallback.
    """
    gene_names = genes["gene"].tolist()
    seed_or_expand = genes["seed_or_expand"].tolist()

    seed_mask = np.array([s == "S" for s in seed_or_expand])
    expand_mask = np.array([s == "E" for s in seed_or_expand])
    seed_idx = np.where(seed_mask)[0]
    expand_idx = np.where(expand_mask)[0]

    full_adj = build_adjacency(gene_names, interactions)
    seed_names = [gene_names[i] for i in seed_idx]
    seed_adj = full_adj[np.ix_(seed_idx, seed_idx)]
    seed_v0 = initial_score_from_adj(seed_adj)
    seed_score = spinner_iteration(seed_adj, seed_v0, max_iter=max_iter, sigma=sigma)

    candidate_names = [gene_names[i] for i in expand_idx]
    raw_pval = expansion_pvalue(
        seed_names,
        candidate_names,
        interactions,
        global_degree,
        total_connected_genes=total_connected_genes,
    )
    q = bh_fdr(raw_pval)
    if q.size and q.min() < expansion_alpha:
        effective_pval = q
    else:
        effective_pval = raw_pval

    keep = effective_pval < expansion_alpha
    kept_candidates = [c for c, k in zip(candidate_names, keep) if k]
    kept_pvals = effective_pval[keep] if keep.any() else effective_pval

    # Iterative expansion — mirrors MATLAB's one-at-a-time addition loop.
    full_nodes = seed_names + kept_candidates
    n_seed = len(seed_names)
    n_full = len(full_nodes)
    full_net = np.zeros((n_full, n_full), dtype=np.float64)
    full_net[:n_seed, :n_seed] = seed_adj

    # Vectorised fill of the seed-to-candidate and candidate-to-candidate
    # blocks. Only edges with at least one candidate endpoint are new; the
    # seed-seed block is already copied from ``seed_adj``.
    full_index = pd.Series(np.arange(n_full, dtype=np.int64), index=full_nodes)
    i1 = interactions["gene1"].map(full_index).to_numpy()
    i2 = interactions["gene2"].map(full_index).to_numpy()
    w = interactions["weight"].to_numpy(dtype=np.float64)
    valid = (~pd.isna(i1)) & (~pd.isna(i2))
    i1v = i1[valid].astype(np.int64)
    i2v = i2[valid].astype(np.int64)
    wv = w[valid]
    touches_candidate = (i1v >= n_seed) | (i2v >= n_seed)
    i1v = i1v[touches_candidate]
    i2v = i2v[touches_candidate]
    wv = wv[touches_candidate]
    full_net[i1v, i2v] = wv
    full_net[i2v, i1v] = wv

    previous_idx = list(range(n_seed))
    in_previous = np.zeros(n_full, dtype=bool)
    in_previous[:n_seed] = True
    previous_nodes = list(seed_names)
    node_origin = ["S"] * n_seed
    previous_score = seed_score

    num_extend = len(kept_candidates)
    iterations_log: list[dict] = []
    neg_inf = np.float64(-np.inf)
    for _ in range(min(num_extend, max_expansions)):
        rank_score = np.zeros(n_full, dtype=np.float64)
        rank_score[previous_idx] = previous_score
        propagated = full_net.T @ rank_score
        # Mask out already-added nodes and argmax in one shot.
        propagated = np.where(in_previous, neg_inf, propagated)
        added = int(np.argmax(propagated))
        if propagated[added] == neg_inf:
            break
        previous_idx.append(added)
        in_previous[added] = True
        previous_nodes.append(full_nodes[added])
        node_origin.append("E")

        sub = full_net[np.ix_(previous_idx, previous_idx)]
        v0 = initial_score_from_adj(sub)
        previous_score = spinner_iteration(sub, v0, max_iter=max_iter, sigma=sigma)
        iterations_log.append(
            {
                "added_gene": full_nodes[added],
                "size": len(previous_idx),
            }
        )

    # Final ranking + null distribution.
    final_net = full_net[np.ix_(previous_idx, previous_idx)]
    final_v0 = initial_score_from_adj(final_net)
    final_score = spinner_iteration(final_net, final_v0, max_iter=max_iter, sigma=sigma)

    n_final = len(previous_nodes)
    if num_random > 0 and n_final > 1 and np.any(final_net != 0):
        resolved = resolve_device(device)
        if progress:
            print(
                f"[winner] random null: n={num_random}, N={n_final}, device={resolved},"
                f" n_jobs={n_jobs}, chunk={chunk}",
                flush=True,
            )
        edge_pool = final_net[np.triu(final_net, k=1) > 0]
        rand_stack = generate_random_networks_batch(
            final_net,
            n=num_random,
            edge_weights=edge_pool,
            seed=random_seed,
            n_jobs=n_jobs,
        )
        v0_stack = initial_score_from_adj_batch(rand_stack)
        random_scores = _batched_spinner(
            rand_stack, v0_stack, max_iter, sigma, resolved, chunk, n_jobs=n_jobs
        ).T  # (N, R)
        ranking_p = ranking_pvalue(final_score, random_scores)
    else:
        ranking_p = np.ones(n_final, dtype=np.float64)

    seed_pval = np.full(n_seed, np.nan, dtype=np.float64)
    exp_pval_final = np.concatenate([seed_pval, kept_pvals])
    if exp_pval_final.size != n_final:
        # If the iterative expansion added fewer than len(kept_candidates)
        # genes, trim from the end.
        exp_pval_final = exp_pval_final[:n_final]

    return WinnerPValueResult(
        gene_names=previous_nodes,
        seed_or_expand=node_origin,
        winner_score=final_score,
        expansion_pval=exp_pval_final,
        ranking_pval=ranking_p,
        iterations=iterations_log,
    )


def devices_report() -> str:
    devs = available_devices()
    return f"available devices: {', '.join(devs)}"
