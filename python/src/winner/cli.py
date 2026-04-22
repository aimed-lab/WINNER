"""Command-line entry points.

Installed scripts:

* ``winner`` — simple mode (``RunWinner.m`` equivalent).
* ``winner-pvalue`` — p-value mode (``RunWinner_withPValue.m`` equivalent).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from . import __version__
from .backend import available_devices
from .io import (
    read_gene_list,
    read_global_degree,
    read_interactions,
    write_winner_result,
    write_winner_result_with_pval,
)
from .pipeline import run_winner, run_winner_with_pvalue


def _common_parser(description: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--gene-list", required=True, type=Path, help="GeneList.txt")
    p.add_argument(
        "--interactions",
        required=True,
        type=Path,
        help="Interaction.txt",
    )
    p.add_argument("-o", "--output", required=True, type=Path, help="Output file path")
    p.add_argument("--max-iter", type=int, default=100, help="Spinner iterations (default 100)")
    p.add_argument("--sigma", type=float, default=0.85, help="Damping factor (default 0.85)")
    p.add_argument("--version", action="version", version=f"winner-net {__version__}")
    return p


def main_simple(argv: list[str] | None = None) -> int:
    parser = _common_parser("WINNER simple mode (no p-values)")
    args = parser.parse_args(argv)

    t0 = time.perf_counter()
    genes = read_gene_list(args.gene_list)
    interactions = read_interactions(args.interactions)
    result = run_winner(genes, interactions, max_iter=args.max_iter, sigma=args.sigma)
    write_winner_result(
        args.output, result.gene_names, result.seed_or_expand, result.winner_score
    )
    dt = time.perf_counter() - t0
    print(f"wrote {args.output} ({len(result.gene_names)} genes) in {dt:.2f}s", file=sys.stderr)
    return 0


def main_pvalue(argv: list[str] | None = None) -> int:
    parser = _common_parser("WINNER with expansion + ranking p-values")
    parser.add_argument(
        "--global-degree",
        required=True,
        type=Path,
        help="AllGeneGloDeg.txt — global degree per gene",
    )
    parser.add_argument(
        "--num-random",
        type=int,
        default=10000,
        help="Number of random networks for the null (default 10000)",
    )
    parser.add_argument(
        "--max-expansions",
        type=int,
        default=50,
        help="Cap on iteratively added expansion genes (default 50)",
    )
    parser.add_argument(
        "--expansion-alpha",
        type=float,
        default=0.05,
        help="Per-candidate expansion p-value cutoff",
    )
    parser.add_argument(
        "--total-connected-genes",
        type=int,
        default=9967,
        help="N for hypergeometric test — depends on the global PPI DB",
    )
    parser.add_argument(
        "--device",
        choices=["auto", "cpu", "cuda", "mps"],
        default="auto",
        help="Backend for the batched random-null spinner",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="joblib workers for random-network generation (-1 = all cores)",
    )
    parser.add_argument(
        "--chunk",
        type=int,
        default=500,
        help="Random-null batch chunk size (memory tuning)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--list-devices",
        action="store_true",
        help="Print detected devices and exit",
    )
    args = parser.parse_args(argv)

    if args.list_devices:
        print(", ".join(available_devices()))
        return 0

    t0 = time.perf_counter()
    genes = read_gene_list(args.gene_list)
    interactions = read_interactions(args.interactions)
    global_degree = read_global_degree(args.global_degree)
    result = run_winner_with_pvalue(
        genes,
        interactions,
        global_degree,
        total_connected_genes=args.total_connected_genes,
        max_expansions=args.max_expansions,
        expansion_alpha=args.expansion_alpha,
        num_random=args.num_random,
        random_seed=args.seed,
        max_iter=args.max_iter,
        sigma=args.sigma,
        device=args.device,
        n_jobs=args.n_jobs,
        chunk=args.chunk,
        progress=True,
    )
    write_winner_result_with_pval(
        args.output,
        result.gene_names,
        result.winner_score,
        result.expansion_pval,
        result.ranking_pval,
    )
    dt = time.perf_counter() - t0
    print(
        f"wrote {args.output} ({len(result.gene_names)} genes, "
        f"{len(result.iterations)} expansions added) in {dt:.2f}s",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main_simple())
