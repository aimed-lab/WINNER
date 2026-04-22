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