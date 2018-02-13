% l0dl:  L0 norm based dictionary learning with acceleration
% D= l0dl(Y, lambda, D0, C0, opts) sovles
%   argmin 1/2||Y-DC||_F^2 + lambda*||C||_0
%     D,C
%   - Input
%       - Y: input data matrix with each column as an observation
%       - lambda: regularization parameter on code sparsity
%       - D0: initial guess of dictionary
%       - C0: initial guess of sparse code
%       - opts: options
%           - n_iter: number of iterations
%           - theta: % maximum step size of C
%           - mu: % maximum step size of D
%   - Output
%       - D: learned dictionary
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference: Chenglong Bao, Hui Ji, Yuhui Quan, Zuowei Shen, 
%L0 norm based dictionary learning by proximal methods with global convergence,
%IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2014
%-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use of this code is free for research purposes only.
%
%Author:  Yuhui Quan
%
%Last Revision: 20-Jun-2014