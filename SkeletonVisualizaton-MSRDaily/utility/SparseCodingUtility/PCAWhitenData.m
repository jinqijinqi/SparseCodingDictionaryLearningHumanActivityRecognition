function [WhitenMatrix WhiteData]=PCAWhitenData(data)
%------------------------------------------------------
%Author: Jin Qi
%Date:   11/24/2014
%Email:  jqichina@hotmail.com
%copyright2014@cnmc
%------------------------------------------------------
%% 
% whiten input data
% input:
%         data: training data
% output:
%         WhitenMatrix: whiten matrix operation
%         WhiteData:     whitened data
%% PCA Whitening
fprintf('\nPCA Whitening\n');

% Remove DC
data = bsxfun(@minus, data, mean(data, 1));

% Remove the "mean" patch
data = bsxfun(@minus, data, mean(data, 2));

% Compute Covariance Matrix and Eigenstuff
cov = data * data' / size(data, 2);
[E,D] = eig(cov);
d = diag(D);

% Sort eigenvalues in descending order
[dsort, idx] = sort(d, 'descend');

% PCA Whitening (and pick top 99% of eigenvalues)
dsum = cumsum(dsort);
dcutoff = find(dsum > 0.99 * dsum(end), 1);
E = E(:, idx(1:dcutoff));
d = d(idx(1:dcutoff));
V = diag(1./sqrt(d+1e-6)) * E';

%% Whiten the data
whiteData = V * data;
WhiteData=whiteData;
WhitenMatrix=V;