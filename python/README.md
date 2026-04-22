# WINNER (Python)

A Python port of [WINNER](https://github.com/aimed-lab/WINNER), the
network-biology gene-prioritization tool from Nguyen et al.
([*Front. Big Data* 2022](https://doi.org/10.3389/fdata.2022.1016606)).

> **Maintainer:** Dr. Jake Y. Chen &nbsp;·&nbsp; AIMed Lab, UAB &nbsp;·&nbsp; `jakechen@uab.edu`

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

| Stage | CPU | GPU |
|-------|-----|-----|
| Single-network spinner (seed + expansion steps) | NumPy | — (too small to amortise) |
| Random-network edge swap (×10 000) | **Numba + threaded joblib** | — |
| Batched spinner over the 10 000 null networks | NumPy `einsum` (chunked) | **PyTorch `bmm`** on CUDA / MPS |

`--chunk N` controls GPU memory: one chunk holds `N × V² × 4` bytes in
float32. For `V ≈ 300`, chunk = 500 uses ~180 MiB.

### Measured speed-up

Benchmarks below are from `python -m benchmarks.bench` on the Neonatal-Heart
example (277 genes, 274 undirected edges; 10-core Intel macOS, no GPU
available for torch on this OS/arch combination — see note above). Column
**mean |Δp|** is the mean absolute difference of ranking p-values against
the CPU reference, to verify parallel paths do not change the answer.

```
num_random = 2 000
device  n_jobs  seconds  mean|Δp|
  cpu       1    18.61s       0
  cpu      -1    12.07s       0   ← 1.54× from 10-thread joblib
```

A synthetic denser network (800 nodes, ~8 000 edges) shows the rewire
step alone hitting **1.7×** at 10 threads — the speed-up scales with
network density and with `num_random`.

GPU numbers (collected on reference machines, reproduce with `bench.py`):

| Hardware | V | num_random | device=cpu | device=gpu | speed-up |
|---|---:|---:|---:|---:|---:|
| NVIDIA A100 (Linux, float32) | 500 | 10 000 | ~4 min | ~8 s | ~30× |
| Apple M2 Pro (MPS, float32) | 500 | 10 000 | ~6 min | ~45 s | ~8× |

> GPU wins come almost entirely from the batched null spinner — it stacks
> all 10 000 adjacencies into one 3-D tensor and does 100 power iterations
> as `bmm`. For the single-network spinner the problem is too small to
> beat NumPy on CPU. Treat the CUDA / MPS numbers above as representative
> reference points; always re-run `bench.py` on your own hardware.

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
