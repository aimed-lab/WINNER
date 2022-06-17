# WINNER
This is the source code for manuscript: Nguyen, T. et al. Gene Prioritization in Network Biology with WINNER.

The code rank the most important genes in a gene lists (seed genes) using their interactions, also expandable to non-seed (expanded) genes. 

To run this code in a simple case, please look at the example in folder NeonatalHeartCaseStudy.

Input files:
GeneList.txt. This file has two columns. The first column is the gene name. The second column indicates whether the genes are seeded or expanded genes.
Interaction.txt. This file has three columns. The first two columns tell which two genes interact, and the third column tells how strong the interaction is. The third column should be normalized between 0 and 1 for better result.

How to run:
After formatting the input files in the same way to folder NeonatalHeartCaseStudy, open Matlab and run file RunWinner.m. The result will be in file winnerResult.txt

The result (file winnerResult.txt):
This file has three column. The first column is the gene name. The second column indicates whether the genes are seeded or expanded genes. The third column, called 'winner score', tells how important (higher means more important) the gene is.
