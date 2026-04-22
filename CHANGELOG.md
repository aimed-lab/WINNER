# Changelog

All notable changes to WINNER (both MATLAB and Python implementations)
are recorded here. Dates are ISO-8601. Versions on this page map 1-to-1
to GitHub tags of the form `vX.Y.Z-py` for the Python port; the MATLAB
implementation predates this changelog and is the reference from the
2022 paper.

## Unreleased

* Repository reorganised: MATLAB files moved under [`matlab/`](matlab/)
  and the Python port stays in [`python/`](python/). Top-level README
  is now an index + decision guide.

## v0.1.1-py — 2026-04-22

[GitHub release](https://github.com/aimed-lab/WINNER/releases/tag/v0.1.1-py)
· [Release notes](python/RELEASE_NOTES.md)

Performance + docs release. No public-API breakage; no change in
ranking p-values at any tolerance.

Headline changes (Python only):

* Four batched-spinner paths with auto-selection per device and network
  density:
  * `spinner_iteration_sparse_batch` (SciPy CSR, threaded joblib)
  * `spinner_iteration_batch` (NumPy `matmul`, BLAS-backed)
  * `spinner_iteration_torch_batch` (PyTorch `bmm`, CPU / CUDA / MPS)
  * `spinner_iteration_torch_sparse_batch` (PyTorch sparse block-diag
    on CUDA)
* Vectorisation cleanup across `io.py`, `pvalue.py`, `pipeline.py` —
  dropped every per-edge Python loop.
* New README sections: "How WINNER works" walkthrough and "Data input
  requirements" column-by-column spec.
* Two new test cases for sparse / dispatcher parity; existing
  MATLAB-reference parity test continues to pass.

Measured impact on the Neonatal-Heart example (V=283, density≈0.4%,
`num_random=2000`, 10-core Intel macOS):

| version | best wall |
|---|---:|
| v0.1.0-py | 15.6 s |
| v0.1.1-py | **11.6 s** |

Isolated spinner phase on a bigger synthetic (V=600, density=1%,
B=1000):

| path | seconds |
|---|---:|
| dense `matmul` | 166.0 |
| **sparse CSR (10 threads)** | **8.0 (20.7× faster)** |

Reference GPU numbers (re-run `python -m benchmarks.bench` on your
hardware):

| hardware | path | B | CPU best | GPU | speed-up |
|---|---|---:|---:|---:|---:|
| NVIDIA A100 | CUDA sparse block-diag | 10 000 | ~4 min | ~6 s | ~40× |
| NVIDIA A100 | CUDA dense `bmm` | 10 000 | ~4 min | ~8 s | ~30× |
| Apple M2 Pro | MPS per-net sparse | 10 000 | ~6 min | ~45 s | ~8× |

## v0.1.0-py — 2026-04-21

[GitHub release](https://github.com/aimed-lab/WINNER/releases/tag/v0.1.0-py)

First Python release. Ports `RunWinner.m` and `RunWinner_withPValue.m`
with numerical parity to the MATLAB reference (`rtol=1e-8` on the
Neonatal-Heart example).

Highlights:

* Installable Python package (`pip install ./python`) exposing the
  `winner` and `winner-pvalue` console scripts plus a Python API.
* Numba-JIT degree-preserving edge swap and threaded joblib over the
  10 000 random-null networks.
* PyTorch-batched dense personalized PageRank on CUDA / Apple MPS via
  `--device auto|cuda|mps|cpu`.
* MIT license; test suite + benchmark harness.

Measured baseline (Neonatal-Heart example, 10-core Intel macOS,
`num_random=2000`):

| n_jobs | seconds |
|---:|---:|
| 1 | 18.6 |
| 10 (threaded) | 12.1 |

## MATLAB reference (pre-2022)

Original MATLAB implementation shipped with the paper. Lives in
[`matlab/`](matlab/). See that folder's README for how to run and for
the worked Neonatal-Heart example.

Citation:

> Nguyen T, Yue Z, Slominski R, Welner R, Zhang J, Chen JY.
> WINNER: A network biology tool for biomolecular characterization and
> prioritization. *Front Big Data.* 2022 Nov 4;5:1016606.
> doi:[10.3389/fdata.2022.1016606](https://doi.org/10.3389/fdata.2022.1016606).
