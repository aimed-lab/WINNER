%% Notes for running Winner with p-value calculation
% - The algorithm and data structure for this option is fairly complex,
% includes multiple components and prebuilt databases. This option should
% be run under the supervision of a trained bioinformatician to make sure
% that the components are installed correctly and the databases are
% up-to-date.

% - The following component files are required:
%      + AllGeneGloDeg.txt. This summarizes the number of gene-gene
% (globally/universally counted) interactions for each gene in HAPPI v.2.0 database
% (http://discovery.informatics.uab.edu/HAPPI/), in 2017. The gene IDs
% are accroding to UniProt (https://www.uniprot.org/). These file
% content should be manually changed/updated, depending on which gene-gene
% database being used as the 'global/universal' and which gene ID system being used.
%       + sym_generate_srand.m. Given an input network, this file creates a
% random network such that the random network has the same nodes, whereas
% each nodes in the random network has the same node degree to the
% corresponding node in the input network. The edges in the random network
% are (completely) different.
%       + spinnerIteration.m. This file produces winner ranking score.

% - The input file format, all are tab-delimited and should have column headers:
%     + file GeneList.txt. This file has two column. The first column is
% the gene identifier, which should be the same to the identification in file
% AllGeneGloDeg.txt. The second column say whether the gene is in the
% 'original' (seed) network (indicated by S) or the expanded gene (indicated by E).
%     + file Interaction.txt. This file has three columns. Columns 1 and 2
% are interactor/interactee genes. These column should use the same gene
% identification to what are used in file GeneList.txt. The third column
% quantifies the interaction strength, normalized between 0 and 1.

% - Output files:
%    + file winnerResult_withPVal.txt. This file has 4 columns: gene
% identification, winner score, 'expansion p-value' and 'ranking
% p-value'. Expansion p-value measures how likely the expanded candidate
% gene (indicated by E in file GeneList.txt) is a good expansion
% (whether p-value < 0.05). Ranking p-value indicates whether the
% ranking for each gene is statistically different from the same gene
% ranking in random networks (randomize preserving node degree)
%    + multiple folders named 'iter <a number from 1 to usually 30>'.
% Expansion and ranking is done by expanding to one gene at an interation.
% Inside each 'iter' folder, there is a 'network file' and 'winner result'
% file being used at the iteration number. The user can manually examine
% them and further decide 'at which iteraction when the expansion is the
% best', then reformat the input files, and get the final answer with the
% RunWinner_simple execution.


%% importdata
geneTable = readtable('GeneList.txt'); % the 'GeneList' file has two columns. The first column is gene name/symbol. The second column is whether the gene is seeded or expanded gene
geneName = table2cell( geneTable(:, 1) );
seedOrExpand = table2cell ( geneTable(:, 2) );
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

%% extract the seed network and winner ranking for seed genes
seedIdx = find( ismember(seedOrExpand, {'S'}) == 1 );
seedGene = geneName(seedIdx);
seedPPI = PPI (seedIdx, seedIdx);

nodeWDeg = sum(seedPPI, 1)';
nodeDeg = sum(sign(seedPPI), 1)';
%initialScore = ones(length(previousNode), 1) / length(previousNode);
initialScore = exp( 2*log(nodeWDeg) - log(nodeDeg));
initialScore(find(isnan(initialScore)==1)) = 0;
[ winnerScore, spinnerIter ] = spinnerIteration( seedPPI, initialScore); % winnerScore stores the winner ranking score

%% ranking, with one-by-one extension
extendIdx = find( ismember(seedOrExpand, {'E'}) == 1 );
extendedCandidate = geneName(extendIdx);

% count the (global) number of interactions for each seed gene
totalDeg = importdata('AllGeneGloDeg.txt');
totalDeg.textdata = strrep(totalDeg.textdata, '_HUMAN', '');
seedGloDeg = zeros(length(seedGene), 1);
for i = 1: length(seedGene)
    [~, index] = ismember(seedGene{i}, totalDeg.textdata);
    if index > 0
        seedGloDeg(i) = totalDeg.data(index);
    end
end

% extendsion p-value
extendPVal = ones(length(extendedCandidate), 1);
N = 9967; % total of gene with connectivity, in UniProtID. This number may be outdated and depends on specific interaction database
n = length(seedGene);

for i = 1 : length(extendPVal)
    index1 = find( ismember(PPITable(:, 1), extendedCandidate(i)) == 1 );
    if index1(1)>0
        k = length( intersect( unique(PPITable(index1, 2)), seedGene) );
        K = totalDeg.data(i);
        extendPVal(i) = 1-hygecdf(k,N,K,n) + hygepdf(k,N,K,n);
    end
end

Q = ones(size(extendPVal));
try
    [FDR, Q] = mafdr(extendPVal);
catch
    Q = ones(size(extendPVal));
end
if min(Q) < 0.05
    extendPVal = Q;
end

if min(extendPVal) < 0.05 % limit the expansion to genes with expansion p-value < 0.05
    extendedCandidate = extendedCandidate(find(extendPVal < 0.05));
    extendPVal = extendPVal(find(extendPVal < 0.05));
end
numExtend = length(extendedCandidate);

%% expansion in each iteration, expanding to one 'best' gene at a time
% iteratively add nodes
FullExNode = [seedGene; extendedCandidate]; % combining the seed genes with the (candidate) extended genes
FullExNetworkAdjust = zeros(length(FullExNode));
FullExNetworkAdjust(1:length(seedGene), 1:length(seedGene)) = seedPPI;
for i = 1 : length(PPITable)
    [~, index1] = ismember(PPITable{i, 1}, extendedCandidate);
    [~, index2] = ismember(PPITable{i, 2}, FullExNode);

    if index1 > 0 && index2 > 0
        FullExNetworkAdjust(index1 + length(seedGene), index2) = PPITable{i, 3};
        FullExNetworkAdjust(index2, index1 + length(seedGene)) = PPITable{i, 3};
    end
end

previousIndex = (1 : length(seedGene))';
previousNode = seedGene;
nodeOrigin = cell(length(seedGene), 1); nodeOrigin(1:length(nodeOrigin)) = {'S'};
previousScore = winnerScore;
for i = 1 : min(numExtend, 50) % maxiumumly adding 50 genes
    try
        rmdir(['iter', num2str(i)]);
    catch
    end
    mkdir(['iter', num2str(i)]);
    rankScore = zeros(length(FullExNode), 1);
    rankScore(previousIndex) = previousScore;
    rankScore = FullExNetworkAdjust' * rankScore;
    [sortRankScore, sortIndex] = sort(rankScore, 'descend');

    candiIndex = setdiff(find(Q<0.05), previousIndex);
    [FullExNode(candiIndex), num2cell(rankScore(candiIndex))];

    for j = 1 : length(sortIndex)
        if ismember(sortIndex(j), previousIndex) == 0
            previousIndex = [previousIndex; sortIndex(j)];
            previousNode = [previousNode; FullExNode(sortIndex(j))];
            nodeOrigin = [nodeOrigin; 'E'];
            break;
        end
    end

    previousNet = FullExNetworkAdjust(previousIndex, previousIndex);
    nodeWDeg = sum(previousNet, 1)';
    nodeDeg = sum(sign(previousNet), 1)';
    initialScore = ones(length(previousNode), 1) / length(previousNode);
    initialScore = exp( 2*log(nodeWDeg) - log(nodeDeg));
    initialScore(find(isnan(initialScore)==1)) = 0;
    [ previousScore, spinnerIter ] = spinnerIteration( previousNet, initialScore);

    [a, b] = find (previousNet ~= 0);
    keepIndex = find(a<b);
    allEdgeScore = zeros(length(a), 1);
    for j = 1 : length(allEdgeScore)
        allEdgeScore(j) = previousNet(a(j), b(j));
    end

    [previousNode, num2cell(previousScore), nodeOrigin];

    xlswrite(['iter', num2str(i), '/winnerUpdate.xlsx'], [previousNode, num2cell(previousScore), nodeOrigin]);
    xlswrite(['iter', num2str(i), '/netUpdate.xlsx'], [previousNode(a), previousNode(b), num2cell(allEdgeScore)]);
end

%% after expansion, rank the final network and compute ranking p-value
finalGeneList = previousNode;
finalNetwork = previousNet;
nodeWDeg = sum(finalNetwork, 1)';
nodeDeg = sum(sign(finalNetwork), 1)';
initialScore = ones(length(finalGeneList), 1) / length(finalGeneList);
initialScore = exp( 2*log(nodeWDeg) - log(nodeDeg));
initialScore(find(isnan(initialScore)==1)) = 0;
[ finalScore, spinnerIter ] = spinnerIteration( finalNetwork, initialScore);

% randomly create 10000 network, preserving node degree for ranking p-value
numExp = 10000;
allRandomScore = zeros(length(finalGeneList), numExp);
netEdgeList = finalNetwork(find(finalNetwork~=0));
for iter = 1:numExp
    randomNet = sym_generate_srand(sign(finalNetwork));
    [c, d] = find(randomNet~=0);
    for i = 1:length(c)
        if c(i) < d(i)
            % pick a random edge weight
            randomWeight = netEdgeList( randi(length(netEdgeList)) );
            randomNet(c(i), d(i)) = randomWeight;
            randomNet(d(i), c(i)) = randomWeight;
        end
    end

    nodeWDeg = sum(randomNet, 1)';
    nodeDeg = sum(sign(randomNet), 1)';
    initialScore = exp( 2*log(nodeWDeg) - log(nodeDeg));
    initialScore(find(isnan(initialScore)==1)) = 0;
    [ allRandomScore(:, iter), ~ ] = spinnerIteration( randomNet, initialScore);
end

rankingPVal = ones(length(finalGeneList), 1);
for i = 1 : length(rankingPVal)
    rankingPVal(i) = length(find( finalScore(i) <= allRandomScore(i, :) )) / size(allRandomScore, 2);
end

expansionPVal = [NaN(length(seedGene), 1); extendPVal];

% write the final result
writetable( table(finalGeneList, finalScore, expansionPVal, rankingPVal), ...
    'winnerResult_withPVal.txt', 'Delimiter', '\t');

