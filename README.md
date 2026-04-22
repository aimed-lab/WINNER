# WINNER
This is the source code for manuscript: Nguyen, T. et al. Gene Prioritization in Network Biology with WINNER.

The code rank the most important genes in a gene lists (seed genes) using their interactions, also expandable to non-seed (expanded) genes.

Two reference implementations are provided:

* **MATLAB** (original) — `spinnerIteration.m` and the `NeonatalHeartCaseStudy/`
  example. See the instructions below.
* **Python** (v0.1.0, April 2026) — a parallel / GPU-enabled port maintained
  by Dr. Jake Y. Chen. Installs with `pip install ./python` and exposes
  `winner` and `winner-pvalue` command-line tools. Supports multi-core CPU
  parallelism and CUDA / Apple-MPS GPU back-ends. Numerical parity with the
  MATLAB reference output on the Neonatal-Heart example is validated by a
  test in the Python package. See [`python/`](python/) and
  [`python/RELEASE_NOTES.md`](python/RELEASE_NOTES.md).

To run the MATLAB code in a simple case, please look at the example in folder NeonatalHeartCaseStudy.

## Citing us:
Please cite the following publication:

Nguyen T, Yue Z, Slominski R, Welner R, Zhang J, Chen JY. WINNER: A network biology tool for biomolecular characterization and prioritization. Front Big Data. 2022 Nov 4;5:1016606. doi: 10.3389/fdata.2022.1016606. PMID: 36407327; PMCID: PMC9672476.


## Input files:
GeneList.txt. This file has two columns. The first column is the gene name. The second column indicates whether the genes are seeded or expanded genes.
Interaction.txt. This file has three columns. The first two columns tell which two genes interact, and the third column tells how strong the interaction is. The third column should be normalized between 0 and 1 for better result.

## How to run (MATLAB):
After formatting the input files in the same way to folder NeonatalHeartCaseStudy, open Matlab and run file RunWinner.m or RunWinner_withPValue.m. Both files use the same input files and format. RunWinner.m is simpler and faster, but it does not show ranking p-value. RunWinner_withPValue.m shows ranking (pr) and expansion (pr) pvalues according to the manuscript above; however, it takes much more time and should be conducted by a trained bioinformatician. Please read the further instructions in file RunWinner_withPValue.m.

## How to run (Python):
```bash
pip install ./python          # or: pip install "./python[all]" for Numba + PyTorch GPU

winner         --gene-list NeonatalHeartCaseStudy/GeneList.txt \
               --interactions NeonatalHeartCaseStudy/Interaction.txt \
               -o winnerResult.txt

winner-pvalue  --gene-list NeonatalHeartCaseStudy/GeneList.txt \
               --interactions NeonatalHeartCaseStudy/Interaction.txt \
               --global-degree NeonatalHeartCaseStudy/AllGeneGloDeg.txt \
               -o winnerResult_withPVal.txt \
               --num-random 10000 --device auto --n-jobs -1
```

The Python version accepts the same input files and produces the same output
format as the MATLAB version. See [`python/README.md`](python/README.md) for
device options, benchmarks, and the Python API.

## The result files:
File winnerResult.txt: The file is the output from RunWinner.m. This file has three column. The first column is the gene name. The second column indicates whether the genes are seeded or expanded genes. The third column, called 'winner score', tells how important (higher means more important) the gene is.
File winnerResult_withPVal.txt: The file is the output from RunWinner.m. This file has four column. The first column is the gene name. The second column, called 'winner score', tells how important (higher means more important) the gene is. The third column is expansion p-value (pe), and the forth column is the ranking p-value (pr) for each gene. Expansion p-value is NaN if the gene is a seeded gene.
