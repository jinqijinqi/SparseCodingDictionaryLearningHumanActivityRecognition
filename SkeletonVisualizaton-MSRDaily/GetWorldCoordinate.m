function GetWorldCoordinate()
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
ProcessedSklDatabase='./Database/MSRAction3DSkeletonReal3D/';
MSRDailyOrigSklDatabase='./MSRDailyUnProcessedSkl/';
%%
mkdir(ProcessedSklDatabase);

DesDir=ProcessedSklDatabase;
SrcDir=MSRDailyOrigSklDatabase;
FileList=dir(fullfile(SrcDir,'*.txt'));
FileList={FileList.name};

for i=1:numel(FileList)
    FileName=FileList{i}
%     if strcmp(FileName,'a13_s06_e01_skeleton.txt');
%         pause;
%     end
    SrcFile=fullfile(SrcDir,FileList{i});
    fh=fopen(SrcFile);
    nFrameJoint=textscan(fh,'%f',2);
    celldisp(nFrameJoint);
    nFrame=nFrameJoint{1};
    nFrame=nFrame(1);
    xyz=[]; 
    for j=1:nFrame
        nRow=cell2mat(textscan(fh,'%d',1))
%         if nRow==0
%             
%             pause;
%         end
        if nRow~=0 % skeleton
        xyzc=textscan(fh,'%f',nRow*4);
        xyzc=reshape(cell2mat(xyzc),4,[]);% each column for coordinate vector (x,y,z,c)
        xyzc=xyzc(:,1:2:40);% odd columns for world coordinates
        xyz=[xyz xyzc(1:3,:)];%   remove confidence level
        end
    end
    fclose(fh);
    DesFile=[DesDir,FileName(1:end-4),'3D.txt'];
    fhDes=fopen(DesFile,'w');
    fprintf(fhDes,'%f %f %f\n',xyz);
    fclose(fhDes);
end

