"""Input / output helpers matching the MATLAB WINNER file formats.

* ``GeneList.txt`` — two tab-delimited columns with header ``Gene`` and
  ``IsSeeded`` (values ``S`` for seed, ``E`` for expansion candidate).
* ``Interaction.txt`` — three tab-delimited columns, ``#node1`` / ``node2`` /
  ``combined_score``. Edge weights are treated as undirected.
* ``AllGeneGloDeg.txt`` — two tab-delimited columns mapping a gene identifier
  (optionally suffixed with ``_HUMAN``) to its global degree.
* ``winnerResult.txt`` — tab-delimited output with three columns.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


def read_gene_list(path: str | Path) -> pd.DataFrame:
    """Read a WINNER ``GeneList.txt`` file, returning a DataFrame with columns
    ``gene`` and ``seed_or_expand`` (values ``S``/``E``)."""
    df = pd.read_csv(path, sep="\t", dtype=str, engine="python")
    df.columns = [c.strip().lstrip("#") for c in df.columns]
    df = df.rename(columns={df.columns[0]: "gene", df.columns[1]: "seed_or_expand"})
    df["gene"] = df["gene"].astype(str).str.strip()
    df["seed_or_expand"] = df["seed_or_expand"].astype(str).str.strip().str.upper()
    return df.reset_index(drop=True)


def read_interactions(path: str | Path) -> pd.DataFrame:
    """Read a WINNER ``Interaction.txt`` file (3 columns)."""
    df = pd.read_csv(path, sep="\t", engine="python")
    df.columns = [c.strip().lstrip("#") for c in df.columns]
    df = df.rename(
        columns={
            df.columns[0]: "gene1",
            df.columns[1]: "gene2",
            df.columns[2]: "weight",
        }
    )
    df["gene1"] = df["gene1"].astype(str).str.strip()
    df["gene2"] = df["gene2"].astype(str).str.strip()
    df["weight"] = pd.to_numeric(df["weight"], errors="coerce").fillna(0.0)
    return df.reset_index(drop=True)


def read_global_degree(path: str | Path) -> dict[str, int]:
    """Return a ``{gene: degree}`` map from ``AllGeneGloDeg.txt``.

    ``_HUMAN`` suffixes are stripped (matches ``RunWinner_withPValue.m``).
    """
    df = pd.read_csv(path, sep="\t", header=None, engine="python")
    out: dict[str, int] = {}
    for _, row in df.iterrows():
        key = str(row[0]).replace("_HUMAN", "").strip()
        try:
            val = int(row[1])
        except (TypeError, ValueError):
            val = 0
        out[key] = val
    return out


def build_adjacency(
    gene_names: Iterable[str],
    interactions: pd.DataFrame,
) -> np.ndarray:
    """Build a symmetric weighted adjacency matrix.

    Only edges where both endpoints appear in ``gene_names`` are kept. If the
    same pair appears multiple times the latest weight wins (matches MATLAB's
    overwriting behaviour in ``RunWinner.m``).
    """
    names = list(gene_names)
    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    adj = np.zeros((n, n), dtype=np.float64)
    g1 = interactions["gene1"].to_numpy()
    g2 = interactions["gene2"].to_numpy()
    w = interactions["weight"].to_numpy(dtype=np.float64)
    for a, b, weight in zip(g1, g2, w):
        i = index.get(a)
        j = index.get(b)
        if i is None or j is None:
            continue
        adj[i, j] = weight
        adj[j, i] = weight
    return adj


def write_winner_result(
    path: str | Path,
    gene_names: Iterable[str],
    seed_or_expand: Iterable[str],
    winner_score: np.ndarray,
) -> None:
    """Write a tab-delimited result file matching ``winnerResult.txt``."""
    df = pd.DataFrame(
        {
            "geneName": list(gene_names),
            "seedOrExpand": list(seed_or_expand),
            "winnerScore": np.asarray(winner_score, dtype=np.float64),
        }
    )
    df.to_csv(path, sep="\t", index=False)


def write_winner_result_with_pval(
    path: str | Path,
    gene_names: Iterable[str],
    winner_score: np.ndarray,
    expansion_pval: np.ndarray,
    ranking_pval: np.ndarray,
) -> None:
    df = pd.DataFrame(
        {
            "finalGeneList": list(gene_names),
            "finalScore": np.asarray(winner_score, dtype=np.float64),
            "expansionPVal": np.asarray(expansion_pval, dtype=np.float64),
            "rankingPVal": np.asarray(ranking_pval, dtype=np.float64),
        }
    )
    df.to_csv(path, sep="\t", index=False, na_rep="NaN")
