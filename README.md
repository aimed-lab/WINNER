# WINNER
This is the source code for manuscript: Nguyen, T. et al. Gene Prioritization in Network Biology with WINNER.

The code rank the most important genes in a gene lists (seed genes) using their interactions, also expandable to non-seed (expanded) genes. 

To run this code in a simple case, please look at the example in folder NeonatalHeartCaseStudy.

## Input files:
GeneList.txt. This file has two columns. The first column is the gene name. The second column indicates whether the genes are seeded or expanded genes.
Interaction.txt. This file has three columns. The first two columns tell which two genes interact, and the third column tells how strong the interaction is. The third column should be normalized between 0 and 1 for better result.

## How to run:
After formatting the input files in the same way to folder NeonatalHeartCaseStudy, open Matlab and run file RunWinner.m or RunWinner_withPValue.m. Both files use the same input files and format. RunWinner.m is simpler and faster, but it does not show ranking p-value. RunWinner_withPValue.m shows ranking (pr) and expansion (pr) pvalues according to the manuscript above; however, it takes much more time and should be conducted by a trained bioinformatician. Please read the further instructions in file RunWinner_withPValue.m.

## The result files:
File winnerResult.txt: The file is the output from RunWinner.m. This file has three column. The first column is the gene name. The second column indicates whether the genes are seeded or expanded genes. The third column, called 'winner score', tells how important (higher means more important) the gene is.
File winnerResult_withPVal.txt: The file is the output from RunWinner.m. This file has four column. The first column is the gene name. The second column, called 'winner score', tells how important (higher means more important) the gene is. The third column is expansion p-value (pe), and the forth column is the ranking p-value (pr) for each gene. Expansion p-value is NaN if the gene is a seeded gene.
