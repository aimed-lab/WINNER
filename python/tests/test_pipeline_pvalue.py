"""Smoke test for the full p-value pipeline.

The p-value pipeline is stochastic (10 k random networks) and so we cannot
compare numerically against the MATLAB reference. Instead, verify the output
has the expected shape / sign / dtype and that seed rows have NaN expansion
p-values.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from winner.io import read_gene_list, read_global_degree, read_interactions
from winner.pipeline import run_winner_with_pvalue


DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize("device", ["cpu"])
def test_pvalue_pipeline_small_null(device):
    genes = read_gene_list(DATA / "GeneList.txt")
    interactions = read_interactions(DATA / "Interaction.txt")
    deg = read_global_degree(DATA / "AllGeneGloDeg.txt")
    result = run_winner_with_pvalue(
        genes,
        interactions,
        deg,
        num_random=50,  # just a smoke test
        device=device,
        n_jobs=1,
        random_seed=7,
        max_expansions=5,
    )
    assert np.isnan(result.expansion_pval[: sum(s == "S" for s in result.seed_or_expand)]).all()
    assert np.all(np.isfinite(result.winner_score))
    assert np.all((result.ranking_pval >= 0) & (result.ranking_pval <= 1))
