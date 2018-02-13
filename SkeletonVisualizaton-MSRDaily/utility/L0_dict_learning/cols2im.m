function [ I ] = cols2im(C, sz_im)
%% cols2im Rearrange matrix columns into an image
%   [ I ] = cols2im(C, im_size) is the reverse operation of im2cols.
%   Input
%       - C: columns 
%       - sz_im: image size
%   Output
%       - I: 2D image matrix
%   See also IM2COL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use of this code is free for research purposes only.
%
%Author:  Yuhui Quan
%
%
%Last Revision: 22-Jun-2014
 %
 %
%%
    sz_step = 1;        %%%SLIDING

    [area_patch, n_patch] = size(C);
    sz_patch = sqrt(area_patch);
 
    I = zeros(sz_im); 
    W = zeros(sz_im); 
    
    i = 1; j = 1;
    for k = 1: 1: n_patch,
        P = reshape(C(:, k), [sz_patch, sz_patch]);
        I(i: i + sz_patch - 1,j: j + sz_patch - 1) = I(i: i + sz_patch - 1,j: j + sz_patch - 1) + P; 
        W(i: i + sz_patch - 1,j: j + sz_patch - 1) = W(i: i + sz_patch - 1, j: j + sz_patch - 1) + 1; 
        if i < sz_im(1) - sz_patch + 1 
            i = i + sz_step; 
        else
            i = 1;
            j = j + sz_step; 
        end;
    end;
    I = I./W;
    
end