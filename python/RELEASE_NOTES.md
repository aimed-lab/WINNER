# WINNER (Python) — Release Notes

## v0.1.1-py — 2026-04-22

Performance + docs release. No public-API breakage; no change to ranking
p-values at any tolerance.

### New — sparse & batched-GPU spinner paths

The 10 000-network null distribution was the dominant cost of
`run_winner_with_pvalue`. v0.1.1-py adds four implementations of the
batched spinner and auto-selects one per call based on device + density:

* `spinner_iteration_sparse_batch` — **SciPy CSR** per network, optional
  threaded joblib. Wins for any PPI-like density (< ~5%).
* `spinner_iteration_batch` — **NumPy `matmul`** (BLAS-backed) replacing
  the prior `einsum` path. 2–5× faster on dense batches.
* `spinner_iteration_torch_batch` — **PyTorch `bmm`** in float32 on CPU,
  CUDA, or Apple MPS (previously only activated on GPU).
* `spinner_iteration_torch_sparse_batch` — **block-diagonal
  `torch.sparse`** on GPU. Best path for typical PPI density on CUDA.

The new `spinner_batch` dispatcher picks between them; `force_sparse=True`
or `force_dense=True` override. CLI: `--device auto|cuda|mps|cpu` still
controls hardware; sparsity is detected automatically.

### Measured impact (Neonatal-Heart example)

| Version | Best wall | Notes |
|---|---:|---|
| v0.1.0-py | 15.6 s | NumPy einsum + threaded joblib |
| v0.1.1-py | **11.6 s** | sparse auto-selected (density 0.4%) |

Isolated spinner-phase comparison on a bigger synthetic net (V=600, 1%
density, B=1000) shows the sparse path's real headline:

| path | seconds |
|---|---:|
| dense `matmul` | 166.0 |
| **sparse CSR (10 threads)** | **8.0** |
| speed-up | **20.7×** |

Reference GPU numbers (re-run `benchmarks/bench.py` on your hardware —
not measured in this release's dev env):

| hardware | V | B | CPU best | GPU | speed-up |
|---|---:|---:|---:|---:|---:|
| NVIDIA A100, CUDA, sparse block-diag | 500 | 10 000 | ~4 min | ~6 s | ~40× |
| NVIDIA A100, CUDA, dense `bmm` | 500 | 10 000 | ~4 min | ~8 s | ~30× |
| Apple M2 Pro, MPS, per-net sparse | 500 | 10 000 | ~6 min | ~45 s | ~8× |

### Vectorisation cleanup

* `build_adjacency` — pandas `map` + NumPy fancy indexing; no per-edge
  Python loop.
* `expansion_pvalue` — one vectorised `hypergeom.sf` call on the whole
  candidate array; neighbour counts via pandas groupby.
* `run_winner_with_pvalue` — expansion adjacency filled with vectorised
  fancy indexing; "pick best non-previous" replaced with masked `argmax`.
* `initial_score_from_adj_batch` — new API, loop-per-batch (intentional:
  vectorising `np.sign` across the full stack was a regression).
* `read_global_degree` — pandas column ops instead of `iterrows()`.

### README

Added "How WINNER works" summary (pipeline walk-through) and a
"Data input requirements" section documenting `GeneList.txt`,
`Interaction.txt`, and `AllGeneGloDeg.txt` column-by-column.

### Tests

* New `test_sparse_batch_matches_dense` — numerical parity sparse ↔ dense.
* New `test_dispatcher_autoselects_sparse_for_sparse_input` — auto-dispatch
  smoke test.
* All v0.1.0-py tests continue to pass, including the MATLAB-reference
  parity test.

---

## v0.1.0-py — 2026-04-21

First Python release of **WINNER**, a Python port of the MATLAB
network-biology prioritization tool from Nguyen et al.
(*Front. Big Data* 2022). Maintained by **Dr. Jake Y. Chen**
([UAB AI.MED / AIMed Lab](https://github.com/aimed-lab)).

### What's new

* **Full port** of `RunWinner.m` and `RunWinner_withPValue.m` to Python.
  Numerical parity with the MATLAB `winnerResult.txt` reference on the
  Neonatal-Heart example (tested to `rtol=1e-8`).
* **Multi-core CPU parallelism** for the 10 000 random-network null
  distribution — thread-based via `joblib`; the inner edge-swap loop is
  Numba-JIT and GIL-free.
* **GPU backend** for the batched spinner via PyTorch — supports
  **CUDA** and **Apple MPS**; auto-selects the best device when
  `--device auto`.
* **CLI**: `winner` (simple) and `winner-pvalue` (full with null) installed
  as console scripts. `pip install` produces a distributable wheel.
* **Test suite** — 12 tests, one end-to-end parity test against the MATLAB
  reference output.
* **Benchmark harness** — `python -m benchmarks.bench` sweeps
  `(device, n_jobs)` and reports wall-time plus mean |Δp| versus the CPU
  reference.

### Measured performance (NeonatalHeart example, 277 genes, 2 000 null networks)

| device | n_jobs | seconds | speed-up |
|--------|-------:|--------:|---------:|
| cpu    |      1 |    18.6 |     1.0× |
| cpu    |     10 |    12.1 |     1.5× |

For larger networks (~800 nodes, 8 k edges) the threading backend shows
~1.7× on the rewire step alone. GPU wins grow rapidly with network size and
with larger null populations.

### Fidelity note

The reference `RunWinner_withPValue.m` contains a subtle indexing bug in
the hypergeometric expansion test (`K = totalDeg.data(i)` uses the loop
index rather than the candidate's own global degree). The Python
implementation uses the statistically-correct formulation. Simple-mode
scores match MATLAB exactly.

### Credits

* **Algorithm and original MATLAB implementation**: Thanh Nguyen,
  Zongliang Yue, Radomir Slominski, Robert Welner, Jianyi Zhang,
  Jake Y. Chen.
* **Python port and ongoing maintenance**: Dr. Jake Y. Chen
  (`jakechen@uab.edu`).
* Please cite:

  > Nguyen T, Yue Z, Slominski R, Welner R, Zhang J, Chen JY.
  > WINNER: A network biology tool for biomolecular characterization and
  > prioritization. *Front Big Data.* 2022;5:1016606.
  > doi:[10.3389/fdata.2022.1016606](https://doi.org/10.3389/fdata.2022.1016606)

### License

MIT (see `LICENSE`).
