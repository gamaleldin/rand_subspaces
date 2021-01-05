%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <subspaces>
% Copyright (C) 2016 Gamaleldin F. Elsayed and John P. Cunningham 
%       (see full notice in README)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demonstration of how to use this code package 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numTimes = 100;
numConds = 2;
numNeus = 60;
allData = randn(numNeus, numConds*numTimes); % create data

D1 = allData(:, 1:100);                      % take some part of the data relevant to some computation
dim = 4;                                     % define dimensionaloty of a subspaces
[w, ~, ~] = svd(D1, 0); 
w = w(:, 1:dim);                             % example of some data aligned subspace (pcs)
Ix = align_ix(w, D1);                         % obtain alignment index from data
numSamples = 10000;                          % define the number of the samples of the empirical random distribution
randIx = sample_rand_subspaces(dim, cov(allData.'), numSamples, 'align_ix', D1);
randIx = vertcat(randIx{:});
pVal = sum(Ix>=randIx)./numSamples;
figure;
hold on
hist(randIx, 0:0.01:1);
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = [0.5 0.5 0.5];
plot(Ix, 0, 'ro', 'markerfacecolor', 'r');
title(['P-value ' num2str(pVal)]);
legend('chance distribution', 'data')
xlabel('alignment index')
ylabel('count')
xlim([0 1])
