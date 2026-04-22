# WINNER (Python) — Release Notes

## v0.1.0 — 2026-04-21

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
