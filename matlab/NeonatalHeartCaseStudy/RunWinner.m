%% importdata
geneTable = readtable('GeneList.txt'); % the 'GeneList' file has two columns. The first column is gene name/symbol. The second column is whether the gene is seeded or expanded gene
geneName = table2cell( geneTable(:, 1) ); seedOrExpand = table2cell ( geneTable(:, 2) );
PPI = zeros(length(geneName)); % this square matrix has interaction among the input gene lists. The order in this matrix should be the same to the order in the input gene list
PPITable = table2cell(readtable('Interaction.txt')); % build the interaction matrix from gene-gene interaction files. Format: <Gene1    Gene2   InteractionScore>
for i = 1 : length(PPITable)
    [~, index1] = ismember(PPITable{i, 1}, geneName);
    [~, index2] = ismember(PPITable{i, 2}, geneName);

    if index1 > 0 && index2 > 0
        PPI(index1, index2) = PPITable{i, 3};
        PPI(index2, index1) = PPITable{i, 3};
    end
end

%% ranking all genes, just use the input gene interaction, without extension
nodeWDeg = sum(PPI, 1)';
nodeDeg = sum(sign(PPI), 1)';
%initialScore = ones(length(previousNode), 1) / length(previousNode);
initialScore = exp( 2*log(nodeWDeg) - log(nodeDeg));
initialScore(find(isnan(initialScore)==1)) = 0;
[ winnerScore, spinnerIter ] = spinnerIteration( PPI, initialScore); % winnerScore stores the winner ranking score
writetable(table(geneName, seedOrExpand, winnerScore), 'winnerResult.txt', 'Delimiter', '\t'); % export the result to file winnerResult.txt

%% this function is the core computation of winner
function [ spinnerScore, spinnerIter ] = spinnerIteration( adj, initialScore, maxIter, sigma )
% Spinner score, by page rank, after iteration. adj stand for interaction
% adjacency matrix, initialScore is the initial score in column vector
% format

if (nargin < 4)
    sigma = 0.85;
    if (nargin < 3)
        maxIter = 100;
    end
end

A = adj;
for i = 1:size(adj, 1)
    if sum(adj(i, :)) > 0
        A(i, :) = adj(i, :) / sum(adj(i, :));
    end
end

spinnerIter = zeros(length(initialScore), maxIter);
spinnerIter(:, 1) = initialScore;

for t = 2 : maxIter
    spinnerIter(:, t) = (1 - sigma) * initialScore + sigma * A' * spinnerIter(:, t-1);
end

spinnerScore = spinnerIter(:, maxIter);
end