function DrawSkl()
%generate mesh
%Author: Jin Qi
%Date:   3/27/2014
%Email:  jqichina@hotmail.com
%copyright2014@gru
%%  
clc;
clear;
global PathAdded;% global variable for adding path
PathAdded=false;

CurDir=pwd;
UtilityPath=[CurDir '/utility']; % path of Utility directory 
addpath(genpath(UtilityPath));
PathAdded=true;
%%
TrueJointDir='./Data/';
ImageDatabase='./Database/MSRAction3DSkeletonReal3D/';
%%
fDir=ImageDatabase
a1=1;a2=2;s1=1;s2=2;e1=1;e2=2;
drawskt(a1,a2,s1,s2,e1,e2,fDir);

