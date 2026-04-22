"""WINNER: Python port of the MATLAB network-biology prioritization tool.

Reference:
    Nguyen T, Yue Z, Slominski R, Welner R, Zhang J, Chen JY.
    WINNER: A network biology tool for biomolecular characterization and prioritization.
    Front Big Data. 2022;5:1016606. doi:10.3389/fdata.2022.1016606
"""

from .core import spinner_iteration, spinner_iteration_batch, initial_score_from_adj
from .io import read_gene_list, read_interactions, build_adjacency, write_winner_result
from .pipeline import run_winner, run_winner_with_pvalue, WinnerResult, WinnerPValueResult
from .backend import get_backend, available_devices

__version__ = "0.1.0"

__all__ = [
    "spinner_iteration",
    "spinner_iteration_batch",
    "initial_score_from_adj",
    "read_gene_list",
    "read_interactions",
    "build_adjacency",
    "write_winner_result",
    "run_winner",
    "run_winner_with_pvalue",
    "WinnerResult",
    "WinnerPValueResult",
    "get_backend",
    "available_devices",
]
