# WINNER — MATLAB implementation

The original MATLAB source for Nguyen et al. (Front. Big Data 2022).

## How to run

After formatting your input files in the same way as
`NeonatalHeartCaseStudy/`, open MATLAB and run either `RunWinner.m` or
`RunWinner_withPValue.m`. Both read the same input files.

* `RunWinner.m` — simpler and faster; emits a WINNER score per gene but
  no p-values.
* `RunWinner_withPValue.m` — also emits expansion and ranking p-values
  but is substantially slower (10 000 random-network null). Should be
  run under the supervision of a trained bioinformatician — see the
  header comments in that file for the extra setup required
  (`AllGeneGloDeg.txt`, correct UniProt-style IDs, etc.).

## Input files

| File | Columns |
|------|---------|
| `GeneList.txt` | `Gene`, `IsSeeded` (`S` = seed, `E` = expansion candidate) |
| `Interaction.txt` | `node1`, `node2`, `combined_score` ∈ [0, 1] |
| `AllGeneGloDeg.txt` | gene id, global degree (p-value mode only) |

## Output files

* `winnerResult.txt` — three columns: gene name, seed-or-expand flag,
  WINNER score (higher = more important).
* `winnerResult_withPVal.txt` — four columns: gene name, WINNER score,
  expansion p-value (NaN for seed rows), ranking p-value.

## Files in this folder

```
matlab/
├── spinnerIteration.m                      # core personalized-PageRank kernel
├── pVal_combineUpdate.xlsx                 # supplementary tables
└── NeonatalHeartCaseStudy/                 # worked example
    ├── RunWinner.m                         # simple mode entry point
    ├── RunWinner_withPValue.m              # p-value mode entry point
    ├── spinnerIteration.m
    ├── sym_generate_srand.m                # degree-preserving network randomiser
    ├── GeneList.txt
    ├── Interaction.txt
    ├── AllGeneGloDeg.txt
    ├── winnerResult.txt                    # reference simple-mode output
    └── winnerResult_withPVal.txt           # reference p-value-mode output
```

See the top-level [README](../README.md) for the MATLAB-vs-Python
decision guide, and the [python/](../python/) directory for the Python
port.
