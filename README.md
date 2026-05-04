# WINNER

WINNER is a network-biology tool for ranking and expanding gene lists
using a personalized-PageRank score over a protein-protein-interaction
(PPI) graph. Given *seed* genes (e.g. GWAS hits, DEGs, curated disease
genes) and a weighted PPI, WINNER returns a "winner score" per gene and
can iteratively add well-connected *expansion* genes that most support
the seed set.

Two reference implementations live in this repository, side-by-side:

| Folder | Language | Status |
|--------|----------|--------|
| [`matlab/`](matlab/) | MATLAB (original) | reference implementation from the 2022 paper; unchanged |
| [`python/`](python/) | Python 3.9+ | parallel + GPU-enabled port; maintained by Dr. Jake Y. Chen |

## Which one should I use?

### Pick **MATLAB** if

* You are reproducing results from the paper exactly and want bit-for-bit
  parity with the published numbers.
* You already have a MATLAB license, live inside MATLAB day-to-day, and
  your networks are small enough that run-time isn't an issue.
* You specifically need `RunWinner_withPValue.m` with its published
  `mafdr` behaviour (the Python port uses `statsmodels`' BH-FDR, which
  is the standard implementation but differs in tie-breaking and
  handling of p-values of 1.0).

### Pick **Python** if

* You want **multi-core CPU and/or GPU scalability** for the 10 000
  random-network null (the expensive part of `withPValue` mode).
* You need to call WINNER from a larger data pipeline, notebook, or
  workflow manager (Snakemake, Nextflow, Airflow, etc.) — the package
  exposes a `winner` / `winner-pvalue` CLI **and** a Python API
  (`run_winner`, `run_winner_with_pvalue`).
* You want to run on a cluster without MATLAB licenses, or in a
  GitHub Action, or in a Docker container.
* You're working with larger networks (V ≳ 500) where the sparse-matrix
  spinner path in v0.1.1-py gives ~20× on the dominant phase.

### Performance quick-take

Measured on the Neonatal-Heart example (V = 283, density ≈ 0.4%,
`num_random = 2000`, 10-core Intel macOS):

| Implementation | Best wall |
|---|---:|
| MATLAB `RunWinner_withPValue.m` | *not measured locally — sequential, typically minutes* |
| Python v0.1.0-py | 15.6 s |
| **Python v0.1.1-py** (sparse auto-selected) | **11.6 s** |

The Python-vs-MATLAB gap widens dramatically for larger / sparser
networks and for the full `num_random = 10 000` setting. See
[`python/README.md`](python/README.md) for details on the sparse and
GPU paths and for a reference GPU-speed-up table.

Ranking p-values match across all paths within floating-point
precision; `tests/test_parity.py` in the Python package verifies this
against the MATLAB reference `winnerResult.txt`.

## Quick start

### MATLAB

```matlab
cd matlab/NeonatalHeartCaseStudy
RunWinner                 % emits winnerResult.txt
RunWinner_withPValue      % emits winnerResult_withPVal.txt (slow)
```

### Python

```bash
pip install ./python                  # core
pip install "./python[all]"           # + Numba + PyTorch (GPU)

winner         --gene-list matlab/NeonatalHeartCaseStudy/GeneList.txt \
               --interactions matlab/NeonatalHeartCaseStudy/Interaction.txt \
               -o winnerResult.txt

winner-pvalue  --gene-list matlab/NeonatalHeartCaseStudy/GeneList.txt \
               --interactions matlab/NeonatalHeartCaseStudy/Interaction.txt \
               --global-degree matlab/NeonatalHeartCaseStudy/AllGeneGloDeg.txt \
               -o winnerResult_withPVal.txt \
               --num-random 10000 --device auto --n-jobs -1
```

Install a specific Python release (tags are preserved):

```bash
pip install "git+https://github.com/aimed-lab/WINNER.git@v0.1.1-py#subdirectory=python"
pip install "git+https://github.com/aimed-lab/WINNER.git@v0.1.0-py#subdirectory=python"
```

## Release history

| Tag | Date | What |
|-----|------|------|
| *unversioned* | 2021–2022 | Original MATLAB implementation (see paper) |
| [`v0.1.0-py`](https://github.com/aimed-lab/WINNER/releases/tag/v0.1.0-py) | 2026-04-21 | First Python release: CPU + GPU parallelism, MATLAB-parity parity test |
| [`v0.1.1-py`](https://github.com/aimed-lab/WINNER/releases/tag/v0.1.1-py) | 2026-04-22 | Sparse + batched-GPU spinner, BLAS-backed dense path, vectorisation pass |

Consolidated notes: [`CHANGELOG.md`](CHANGELOG.md). Per-release notes
also live in [`python/RELEASE_NOTES.md`](python/RELEASE_NOTES.md) and on
the [GitHub Releases page](https://github.com/aimed-lab/WINNER/releases).

## Data input requirements (both implementations)

All inputs are **tab-delimited** text with a header row.

### `GeneList.txt`

| column | meaning |
|---|---|
| `Gene` | gene identifier (symbol or UniProt; must match the interaction and global-degree files) |
| `IsSeeded` | `S` = seed gene; `E` = expansion candidate to be scored |

### `Interaction.txt`

| column | meaning |
|---|---|
| `node1` | gene identifier |
| `node2` | gene identifier |
| `combined_score` | edge weight, normalised to `[0, 1]` for best results |

Undirected — listing an edge once is enough. If both directions appear,
the later-seen weight wins.

### `AllGeneGloDeg.txt` (p-value mode only)

| column | meaning |
|---|---|
| gene id | same namespace; a `_HUMAN` suffix is stripped |
| global degree | number of interactions for the gene in the *whole* PPI database (not just your subnet). If you change PPI databases, regenerate this file. |

Example files live in
[`matlab/NeonatalHeartCaseStudy/`](matlab/NeonatalHeartCaseStudy/)
and, copied as test fixtures, in
[`python/tests/data/`](python/tests/data).

## Citation

Please cite the original paper regardless of which implementation you
use:

> Nguyen T, Yue Z, Slominski R, Welner R, Zhang J, Chen JY.
> WINNER: A network biology tool for biomolecular characterization and
> prioritization. *Front Big Data.* 2022 Nov 4;5:1016606.
> doi:[10.3389/fdata.2022.1016606](https://doi.org/10.3389/fdata.2022.1016606).
> PMID: 36407327; PMCID: PMC9672476.

## License

The MATLAB reference retains its existing notices and original authorship in
[`LICENSE`](LICENSE). The Python port is distributed under the non-commercial
research and education license in [`python/LICENSE`](python/LICENSE):
commercial use requires a separate written license granted by Jake Y. Chen,
AIMed Lab, UAB, or another authorized copyright holder.
