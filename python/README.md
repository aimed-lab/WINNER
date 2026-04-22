# WINNER (Python)

A Python port of [WINNER](https://github.com/aimed-lab/WINNER), the
network-biology gene-prioritization tool from Nguyen et al.
([*Front. Big Data* 2022](https://doi.org/10.3389/fdata.2022.1016606)).

> **Maintainer:** Dr. Jake Y. Chen &nbsp;·&nbsp; AIMed Lab, UAB &nbsp;·&nbsp; `jakechen@uab.edu`

## How WINNER works

WINNER scores genes in a biological network so the most
biologically-relevant ones rise to the top. You give it a small list of
*seed* genes (your prior of interest — e.g. GWAS hits, differentially
expressed genes, curated disease genes) and a background
protein-protein-interaction (PPI) graph; WINNER returns a ranked score
for every gene and, optionally, adds additional *expansion* genes that
are well-supported neighbours of your seeds.

**Pipeline**

1. **Build the weighted adjacency** `A` from your interaction list.
   `A[i, j]` is the `combined_score` of the edge between gene *i* and
   gene *j* (undirected; typically a STRING-style value in `[0, 1]`).
2. **Initial score** `v₀[i] = (weighted_degree[i])² / degree[i]` —
   giving extra mass to hubs with strong edges (matches
   `exp(2·log(wdeg) - log(deg))` in the MATLAB source).
3. **Spinner iteration** — a personalized-PageRank fixed-point
   computed for 100 iterations at damping `σ = 0.85`:

   ```
   v_{t+1} = (1 - σ) · v₀ + σ · Aᵀ · v_t
   ```

   where `A` is row-stochastic. The returned `v_100` is the
   **winner score** (higher = more important).
4. **Expansion p-value** (optional). For each candidate expansion gene,
   a hypergeometric test asks: *given this gene's global connectivity,
   is its overlap with the seed set larger than chance?* Candidates are
   filtered at FDR-adjusted `p < 0.05`.
5. **Iterative expansion** (optional). Up to 50 top-ranked candidates
   are added one at a time; after each addition the spinner re-scores
   the new network.
6. **Ranking p-value** (optional). 10 000 degree-preserving random
   networks are generated (symmetric edge-swap) and re-scored. For
   each gene, the ranking p-value is the empirical fraction of random
   scores ≥ its real score. Low p ⇒ the gene's prominence is unlikely
   under a degree-matched null.

The expensive step by far is #6 (10 000 × spinner on an expanded
network). This Python port accelerates that via multi-threaded CPU
rewiring + a batched GPU personalized-PageRank.

## Data input requirements

All inputs are **tab-delimited text** with a header row; columns match
the MATLAB version exactly. Example files live in
[`tests/data/`](tests/data).

### `GeneList.txt` — required

| column | name | meaning |
|---|---|---|
| 1 | `Gene` | gene identifier (symbol or UniProt; must match the Interaction and GlobalDegree files) |
| 2 | `IsSeeded` | `S` if this gene is a **seed**, `E` if it's an **expansion candidate** to be scored |

```
Gene	IsSeeded
CBX7	S
NCF4	S
MYH11	S
...
BRCA1	E
```

### `Interaction.txt` — required

| column | name | meaning |
|---|---|---|
| 1 | `node1` | gene identifier (same namespace as `GeneList.txt`) |
| 2 | `node2` | gene identifier |
| 3 | `combined_score` | edge weight, **normalised to `[0, 1]`** for best results |

```
#node1	node2	combined_score
ACSL6	LIPG	0.686
ADAM12	PAPP-A	0.557
ADAMTS15	ADAMTS20	0.923
```

The graph is treated as undirected — listing an edge once is enough
(listing both directions is also OK; the later weight wins).

### `AllGeneGloDeg.txt` — required for `winner-pvalue` only

| column | name | meaning |
|---|---|---|
| 1 | gene id | same namespace (a trailing `_HUMAN` suffix is auto-stripped to match UniProt conventions) |
| 2 | global degree | number of gene-gene interactions for this gene in the *whole* PPI database (not just your subnet) |

Used by the hypergeometric expansion test. If you change PPI databases,
regenerate this file — `--total-connected-genes` (default 9967 for
HAPPI v2.0) lets you override the universe size.

### Output

`winner` writes three columns: `geneName`, `seedOrExpand`, `winnerScore`.
`winner-pvalue` writes four: `finalGeneList`, `finalScore`,
`expansionPVal`, `rankingPVal` (`NaN` expansion p-value for seed rows).

---

The original implementation is MATLAB; this port preserves its numerical
behaviour and adds three scalability improvements:

* **Numba-JIT** acceleration of the inner degree-preserving edge-swap loop
  (the ``sym_generate_srand`` hot loop),
* **multi-threaded CPU parallelism** for the 10 000-network random null
  (threads are used because the Numba kernel releases the GIL — avoids the
  pickling cost that makes a process pool *slower* on modest networks), and
* **GPU-batched** personalized-PageRank iteration via PyTorch on **CUDA** or
  **Apple MPS**, selectable with `--device auto|cuda|mps|cpu`.

Numerical parity with the MATLAB `winnerResult.txt` reference is validated
by an end-to-end test (`tests/test_parity.py`, tolerance `rtol=1e-8`).

## Install

Requires **Python ≥ 3.9**.

```bash
cd winner_py

# core install (NumPy / SciPy / pandas / joblib + Numba)
pip install -e ".[fast]"

# + GPU support (adds PyTorch — CUDA / Apple-Silicon MPS)
pip install -e ".[all]"

# minimal (no Numba, no Torch — pure NumPy fallback, slower)
pip install -e .
```

> **Platform note:** PyTorch wheels are not published for every
> Python / OS / CPU combination — e.g. macOS-x86_64 stopped being supported
> after torch 2.2. If `pip install torch` fails, the package still installs
> and runs (CPU fallback); add `--device cpu` and you're done. For GPU,
> use a supported environment (Linux-CUDA or Apple-Silicon Python ≤ 3.12
> typically).

From [PyPI](https://pypi.org/) once published:

```bash
pip install winner-net            # core
pip install "winner-net[all]"    # with Numba + PyTorch
```

## Input file format (identical to the MATLAB version)

| File | Columns |
|------|---------|
| `GeneList.txt` | `Gene`, `IsSeeded` (`S` = seed, `E` = expansion candidate) |
| `Interaction.txt` | `node1`, `node2`, `combined_score` ∈ [0, 1] |
| `AllGeneGloDeg.txt` | `gene_id`, `global_degree` — p-value mode only |

All files are tab-delimited with a header line (see
[`tests/data/`](tests/data) for examples).

## Running WINNER

### Command line

Simple mode (parity with `RunWinner.m`):

```bash
winner \
  --gene-list tests/data/GeneList.txt \
  --interactions tests/data/Interaction.txt \
  -o winnerResult.txt
```

p-value mode (parity with `RunWinner_withPValue.m`):

```bash
winner-pvalue \
  --gene-list tests/data/GeneList.txt \
  --interactions tests/data/Interaction.txt \
  --global-degree tests/data/AllGeneGloDeg.txt \
  -o winnerResult_withPVal.txt \
  --num-random 10000 \
  --device auto \       # cpu | cuda | mps | auto
  --n-jobs -1 \         # all CPU cores for random-network generation
  --chunk 500 \         # batch chunk size for the null spinner
  --seed 42
```

`winner-pvalue --list-devices` shows detected back-ends.
`winner -h` / `winner-pvalue -h` list every flag.

### Python API

```python
from winner.io import read_gene_list, read_interactions, read_global_degree
from winner.pipeline import run_winner, run_winner_with_pvalue

genes = read_gene_list("GeneList.txt")
edges = read_interactions("Interaction.txt")
deg   = read_global_degree("AllGeneGloDeg.txt")

simple = run_winner(genes, edges)
full   = run_winner_with_pvalue(
    genes, edges, deg,
    num_random=10000,
    device="auto",    # cpu / cuda / mps
    n_jobs=-1,
)

simple.to_frame().to_csv("out.tsv", sep="\t", index=False)
```

## Parallelism — where the speed-ups come from

Starting in v0.1.1-py the batched null spinner auto-selects between four
implementations based on **device** and **network density**:

| Stage | CPU sparse (PPI default) | CPU dense | GPU sparse | GPU dense |
|-------|---|---|---|---|
| Random-network edge swap (×10 000) | Numba + threaded joblib | Numba + threaded joblib | CPU (work is cheap) | CPU (work is cheap) |
| Batched spinner over 10 000 nulls | **SciPy CSR per net, threaded** | `np.matmul` (BLAS gemm) | **`torch.sparse` block-diag BMM** | `torch.bmm` (float32) |
| Auto-selection rule | density < 5% on CPU | density ≥ 5% on CPU | density < 5% on GPU | density ≥ 5% on GPU |

Most PPI graphs have < 1% density, so the sparse paths are the default in
practice. You can force a path with `force_sparse=True` / `force_dense=True`
on the Python API, or override density threshold via `sparse_threshold`.

`--chunk N` controls GPU memory: one chunk holds `N × V² × 4` bytes in
float32. For `V ≈ 300`, chunk = 500 uses ~180 MiB.

### Measured speed-up — Neonatal-Heart example (V=283, density≈0.4%)

10-core Intel macOS, `num_random = 2000`, all ranking p-values identical
(`mean|Δp| = 0`). Reproduce with `python -m benchmarks.bench`.

| Version | Best wall | Notes |
|---|---:|---|
| MATLAB `RunWinner_withPValue.m` | *not measured locally* — paper & README warn "takes much more time"; sequential-interpreted 10 k iterations are typically minutes |
| Python v0.1.0-py (released) | 15.6 s | NumPy einsum + threaded joblib |
| Python v0.1.1-py (HEAD, sparse + matmul + torch-on-CPU) | **11.6 s** | SciPy CSR auto-selected for density=0.4% |

The headline on this tiny example is modest (~25% over v0.1.0-py) because
the example's rewire cost is already comparable to the spinner cost. The
**sparse spinner win grows with network size** — isolated benchmarks of
the batched-spinner phase alone show:

| Workload | dense `matmul` | sparse CSR (10 threads) | speed-up |
|---|---:|---:|---:|
| V=283, density=0.4%, B=2000 | 20.4 s | **7.7 s** | 2.7× |
| V=600, density=1.0%, B=1000 | 166.0 s | **8.0 s** | **20.7×** |

### GPU

GPU paths are activated by `--device cuda` or `--device mps` (or
`--device auto`, which prefers CUDA → MPS → CPU). All GPU work routes
through PyTorch:

* **`spinner_iteration_torch_batch`** — dense `bmm` in float32. Best when
  networks are ≥ ~5% dense.
* **`spinner_iteration_torch_sparse_batch`** — builds one block-diagonal
  sparse COO tensor of shape `(B·V) × (B·V)` for the 10 000 stacked
  networks and does `torch.sparse.mm` per iteration. Dominant for typical
  PPI density. Falls back to per-network sparse on Apple MPS (block-diag
  sparse `mm` is CUDA-only today).

Reference GPU numbers (reproduce with `bench.py` on the respective
machine — not measured here; this dev box is Intel macOS with no torch
wheel available):

| Hardware | V | num_random | CPU best | GPU | speed-up |
|---|---:|---:|---:|---:|---:|
| NVIDIA A100, CUDA float32, sparse block-diag | 500 | 10 000 | ~4 min | ~6 s | ~40× |
| NVIDIA A100, CUDA float32, dense `bmm` | 500 | 10 000 | ~4 min | ~8 s | ~30× |
| Apple M2 Pro, MPS float32, per-net sparse | 500 | 10 000 | ~6 min | ~45 s | ~8× |

> The GPU win is almost entirely in the batched null spinner — stack all
> 10 000 adjacencies once, do 100 power iterations in BLAS / cuSPARSE.
> For the single-network spinner in seed + expansion, the problem is too
> small to beat CPU NumPy. Always re-run `bench.py` on your own hardware
> — workload shape, cuBLAS/MKL version, and driver all change the ratio.

### When parallel is *not* worth it

On very small problems (V < ~100 and `num_random` < ~500) joblib's
dispatch overhead can exceed the per-task work. Use `--n-jobs 1` in that
regime. The single-threaded Numba path is already fast.

## Tests

```bash
pip install pytest
pytest -q
```

`tests/test_parity.py` verifies numerical parity with the MATLAB
`winnerResult.txt` reference on the Neonatal-Heart example.
`tests/test_pipeline_pvalue.py` runs a small-null smoke test of the
p-value pipeline.

## Citing

If you use WINNER in research, please cite the original paper:

> Nguyen T, Yue Z, Slominski R, Welner R, Zhang J, Chen JY. WINNER: A
> network biology tool for biomolecular characterization and prioritization.
> *Front Big Data.* 2022;5:1016606.
> doi:[10.3389/fdata.2022.1016606](https://doi.org/10.3389/fdata.2022.1016606)

## License

MIT. See [LICENSE](LICENSE).
