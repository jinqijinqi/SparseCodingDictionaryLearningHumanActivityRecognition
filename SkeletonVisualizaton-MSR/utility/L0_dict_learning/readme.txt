README

Version 1.1, 04-July-2014.

This package contains the MATLAB implementation of L0 norm based dictionary learning for sparse coding

==========================================
Use of this code is free for research purposes only.

================================================================================================
Reference:

Chenglong Bao, Hui Ji, Yuhui Quan, Zuowei Shen, L0 norm based dictionary learning by proximal methods with global convergence, CVPR, 2014

=================================================================
This implementation has been tested with MATLAB 2011b on a PC with 64-bit Ubuntu Linux 12.04 LTS.

==================================================================
Installation and additional toolbox requirement

1. Unpack the contents of the compressed file to a new directory.

2. (optional) to running the demo of image denoising, the following two toolboxes are required

	a. the omp2.m from the OMP toolbox  (http://www.cs.technion.ac.il/~ronrubin/software.html) 
		(the compiled package for linux is included for convenience, you may need to re-compile it for solving possible compatibility issue.)
	
	b. The image processing toolbox of MATLAB
. 

======================
Demo

See demo.m for a demonstration of using the learned dictionary for image denoising.

======================================================================
Main Routine


l0dl  learns a dictionary from training samples by the accelerated proximal method


D= l0dl(Y, lambda, D0, C0, opts) solves
 $argmin 1/2||Y-DC||_F^2 + lambda*||C||_0$
   D,C
  
 - Input
      - Y: input data matrix with each column as an observation
      - lambda: regularization parameter on code sparsity
      - D0: initial guess of dictionary
      - C0: initial guess of sparse code
      - opts: options
          - n_iter: number of iterations
          - theta: algorithm parameter, adaptive if not set
          - mu: algorithm parameter,  adaptive if not set
         
 - Output
      - D: learned dictionary
    


