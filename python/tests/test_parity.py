"""End-to-end parity test against the MATLAB reference output.

We run the simple WINNER pipeline on the NeonatalHeart example and compare to
``winnerResult.txt``. The spinner is deterministic, so scores should match to
floating-point tolerance.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from winner.io import read_gene_list, read_interactions
from winner.pipeline import run_winner


DATA = Path(__file__).parent / "data"


def test_run_winner_matches_matlab_reference():
    genes = read_gene_list(DATA / "GeneList.txt")
    interactions = read_interactions(DATA / "Interaction.txt")
    result = run_winner(genes, interactions)

    ref = pd.read_csv(DATA / "winnerResult.txt", sep="\t")
    ref.columns = [c.strip() for c in ref.columns]

    # Align by gene name.
    got = pd.DataFrame(
        {
            "geneName": result.gene_names,
            "winnerScore": result.winner_score,
        }
    )
    merged = got.merge(ref, on="geneName", suffixes=("_py", "_m"))
    assert len(merged) == len(ref), "every reference gene must appear in Python output"

    py = merged["winnerScore_py"].to_numpy()
    m = merged["winnerScore_m"].to_numpy()
    # MATLAB prints ~14 significant digits; allow loose tolerance.
    np.testing.assert_allclose(py, m, rtol=1e-8, atol=1e-10)


def test_run_winner_scores_are_finite_and_nonneg():
    genes = read_gene_list(DATA / "GeneList.txt")
    interactions = read_interactions(DATA / "Interaction.txt")
    result = run_winner(genes, interactions)
    assert np.all(np.isfinite(result.winner_score))
    assert np.all(result.winner_score >= -1e-12)
