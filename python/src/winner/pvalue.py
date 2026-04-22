"""Expansion and ranking p-value computations.

Ports the hypergeometric-test + BH-FDR block and the empirical ranking-p-value
block from ``RunWinner_withPValue.m``.

Note on MATLAB fidelity: the reference ``RunWinner_withPValue.m`` contains a
subtle indexing bug in the expansion test (``K = totalDeg.data(i)`` uses the
loop index ``i`` instead of the candidate's own global degree). We reproduce
the statistically-correct formulation here — the candidate's own global degree
from ``AllGeneGloDeg.txt`` is used for ``K``.
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def expansion_pvalue(
    seed_genes: Sequence[str],
    candidate_genes: Sequence[str],
    interactions: pd.DataFrame,
    global_degree: Mapping[str, int],
    total_connected_genes: int = 9967,
) -> np.ndarray:
    """Hypergeometric expansion p-values for each candidate gene.

    For candidate *c*:
        N = ``total_connected_genes``
        K = global degree of *c*
        n = number of seed genes
        k = number of distinct seed genes *c* interacts with

    Returns ``P(X >= k)`` for each candidate.
    """
    seed_set = set(seed_genes)
    n = len(seed_set)
    candidates = list(candidate_genes)
    if not candidates:
        return np.zeros(0, dtype=np.float64)

    neighbours_of: dict[str, set[str]] = {c: set() for c in candidates}
    cand_set = set(candidates)
    g1 = interactions["gene1"].to_numpy()
    g2 = interactions["gene2"].to_numpy()
    for a, b in zip(g1, g2):
        if a in cand_set:
            neighbours_of[a].add(b)
        if b in cand_set:
            neighbours_of[b].add(a)

    pvals = np.ones(len(candidates), dtype=np.float64)
    for i, c in enumerate(candidates):
        K = int(global_degree.get(c, 0))
        if K <= 0:
            continue
        hits = len(neighbours_of[c] & seed_set)
        pvals[i] = float(hypergeom.sf(hits - 1, total_connected_genes, K, n))
    return pvals


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values (matches MATLAB ``mafdr``)."""
    from statsmodels.stats.multitest import multipletests

    pvals = np.asarray(pvals, dtype=np.float64)
    if pvals.size == 0:
        return pvals.copy()
    _, q, _, _ = multipletests(pvals, method="fdr_bh")
    return q


def ranking_pvalue(
    final_score: np.ndarray,
    random_scores: np.ndarray,
) -> np.ndarray:
    """Empirical one-sided ranking p-value per gene.

    ``P_i = #{r: finalScore_i <= random_scores[i, r]} / R``, matching
    ``length(find( finalScore(i) <= allRandomScore(i, :) )) / size(...)``.
    """
    final_score = np.asarray(final_score, dtype=np.float64)
    rand = np.asarray(random_scores, dtype=np.float64)
    if rand.ndim != 2 or rand.shape[0] != final_score.shape[0]:
        raise ValueError("random_scores must be (N, R) matching final_score")
    hits = (rand >= final_score[:, None]).sum(axis=1)
    return hits / rand.shape[1]
