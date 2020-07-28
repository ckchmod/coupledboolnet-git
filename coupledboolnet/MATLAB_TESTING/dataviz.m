clear all ;clc;
load('test.mat')
s = statetransition(:,1)+1;
t = statetransition(:,2)+1;

%s = [1, 2,3]
%t = [2, 1, 3]
s = cast(s,'double')
t = cast(t,'double')

% Return a Steady State Distributions
numberCells = 1;
numberGenes = 11;
numberT = size(allstates,1)
decMatrix = zeros(1, numberT);

% Convert states to decimal representation
for tt = 1:numberT
    tempDec = 0;
    for jj = 1:numberGenes
        tempDec = tempDec + allstates(tt,jj)*2^(numberGenes-jj);
    end
    decMatrix(1,tt) = tempDec;
end

% Return the Steady State Probability Distribution
counts = zeros(numberCells, 2^numberGenes);
for i = 1:length(decMatrix(:,1))
    for j = 1:length(decMatrix(1,:))
        counts(i, decMatrix(i,j)+1) = counts(i, decMatrix(i,j)+1) + 1;
    end
end
rhs = counts/(sum(counts(1,:)));
figure(1)
set(gca, 'YScale', 'log')
s1 = bar(rhs)
str=sprintf('Myeloid Progenitors Steady-State Distribution: Internal Random Noise p = 0.1, GATA-2 Noise q = .99');
title(str); xlim([0 2048]);

G = digraph(s, t);
figure(2)
plot(G,'Layout','force');
title('State Transition of Hierarchical Differentiation of Myeloid Progenitors');