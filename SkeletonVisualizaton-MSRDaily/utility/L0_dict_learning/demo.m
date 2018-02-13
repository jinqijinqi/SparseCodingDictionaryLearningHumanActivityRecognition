%% This demo reads an image, adds random white noise and denoises it using the dictionary learnt by L0DL_AC and OMP.
disp('====================================================================');
disp('This demo show how to learn an over-complete dictionary from a noisy image,')
disp('and then denoising it using the OMP method on the learned dictionary in one pass')
disp('====================================================================');
%%
clear all; close all; warning('off');
addpath('omp_box');
randn('seed',0); rand('seed',0);
disp('initialization...');

%% read clear image and synthesize the noisy one
im_name = 'lena512.png';
std_noise = 25;     % standard deviation of Gaussian noise
im_gt = double(imread(im_name));    % ground truth image
im_n = im_gt + randn(size(im_gt)) * std_noise;  % noisy image with Gaussian noise
disp(['noisy image =  ', im_name, ' + Gaussian noise w/ s.t.d. = ', num2str(std_noise)]);

%% set parameters %%
sz_patch = 8;   % patch size,  %%%%NOTE THAT ONLY SQUARE PATCH IS SUPPORTED IN THIS VERSION!
n_atom = 256;       % number of atoms in dictionary, i.e., dictioanry size
n_sample = 40000;   % number of training samples
n_iter = 30;    % maxiaum number of iterations in dictionary learning
lambda = 6500;  % model parameter cotrolling sparisty regularization
epsilon = sqrt(sz_patch^2) * std_noise * 1.1;   % target error for omp

%% denoising
% generate training data
A = im2cols(im_n,sz_patch);
idx = randperm(size(A,2));
Y = A(:,idx(1:n_sample));
% A = A - repmat(mean(A,1),[size(A,1) 1]);  %%%ZERO-MEAN ALIGNMENT.NOT BETTER!!

% Setting parameters  
opts_dl.theta = 1;  % maximum step size of C
opts_dl.mu = 1e-3;  % maximum step size of D
opts_dl.n_iter = n_iter;

% Setting the initialization for dictionary learning
D0 = dct2dict(sz_patch,n_atom);     % initial dictionary
C0 = omp2(D0'*Y,sum(Y.^2,1),D0'*D0,epsilon);    % initial sparse code

% dictionary learning
disp('====================================================================');
disp('learning dictionary...');
%tic
D = l0dl(Y, lambda, D0, C0, opts_dl);
%toc

% generating the final sparse code using the  OMP method for denoising
disp('====================================================================');
disp('denoising...');
Y = A;
M = repmat(mean(Y,1),[size(Y,1) 1]);
Y = Y - M;
C = omp2(D'*Y,sum(Y.^2,1),D'*D,epsilon,'maxatoms',floor(sz_patch^2/2),'checkdict','off'); 

% generating denoised image
im_r = cols2im(D*C+M, size(im_gt));

%% output results
disp('====================================================================');
disp(['The dictonary is stored in variable D.']);
figure; subplot(1,2,1); 
imshow(uint8(im_n));
title('noisy image');
subplot(1,2,2);
imshow(uint8(im_r));
title(['de-noised result with psnr value ', num2str(10 * log10(255^2 /mean((im_r(:)-im_gt(:)).^2)))]);
figure; dictshow(D); title('learned dictionary');
