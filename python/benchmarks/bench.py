"""Benchmark WINNER p-value pipeline across devices and parallel configs.

Runs the NeonatalHeart example with varying ``num_random``, ``device``, and
``n_jobs`` and prints a comparison table. Correctness is sanity-checked by
computing the mean absolute difference between each configuration's ranking
p-values and a reference CPU/single-job run on the same seed.

Usage:
    python -m benchmarks.bench                    # default sweep
    python -m benchmarks.bench --num-random 2000  # override null size
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from winner.backend import available_devices
from winner.io import read_gene_list, read_global_degree, read_interactions
from winner.pipeline import run_winner_with_pvalue


DATA = Path(__file__).resolve().parent.parent / "tests" / "data"


def run_one(genes, interactions, deg, *, num_random, device, n_jobs, seed=12345):
    t0 = time.perf_counter()
    res = run_winner_with_pvalue(
        genes,
        interactions,
        deg,
        num_random=num_random,
        random_seed=seed,
        device=device,
        n_jobs=n_jobs,
        progress=False,
    )
    return time.perf_counter() - t0, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--num-random", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    genes = read_gene_list(DATA / "GeneList.txt")
    interactions = read_interactions(DATA / "Interaction.txt")
    deg = read_global_degree(DATA / "AllGeneGloDeg.txt")

    print(f"available devices: {available_devices()}")
    print(f"num_random = {args.num_random}, seed = {args.seed}")
    print()

    configs = [
        ("cpu", 1),
        ("cpu", -1),
    ]
    for dev in ("cuda", "mps"):
        if dev in available_devices():
            configs.append((dev, -1))

    ref_p = None
    header = f"{'device':>6} {'n_jobs':>7} {'seconds':>10} {'mean|Δp|':>12}"
    print(header)
    print("-" * len(header))
    for device, n_jobs in configs:
        dt, res = run_one(
            genes, interactions, deg,
            num_random=args.num_random, device=device, n_jobs=n_jobs, seed=args.seed,
        )
        if ref_p is None:
            ref_p = res.ranking_pval
            diff = 0.0
        else:
            diff = float(np.mean(np.abs(res.ranking_pval - ref_p)))
        print(f"{device:>6} {n_jobs:>7} {dt:>10.3f} {diff:>12.2e}")


if __name__ == "__main__":
    main()
