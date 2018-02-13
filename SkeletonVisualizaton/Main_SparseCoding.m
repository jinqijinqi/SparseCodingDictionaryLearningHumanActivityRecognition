function Main_SparseCoding_Dictionary()
%generate mesh
%Author: Jin Qi
%Date:   1/30/2014
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
JointDatabase='./Result/Joint/';
JointFileName='JointPoint.mat';
% JointDatabase='./Result/JointAngle/';
% JointFileName='_JointPoint.mat';

ResultDir='./Result/';
ResultBasisDir='./Result/Basis/';
ResultIndex='./Result/Index/';
mkdir(ResultDir);
mkdir(ResultBasisDir);

SubjectList=dir(JointDatabase);
SubjectList=SubjectList(3:end); % remove . or .. dir
% load(fullfile(ResultIndex,'File2ClassMap.mat'));
% load(fullfile(ResultIndex,'SubjectFile2ClassMap.mat'));
load(fullfile(ResultIndex,'IndexTable.mat'))
SubjectList=cell2mat([IndexTable(:,2)]);
ClassList=cell2mat([IndexTable(:,3)]);
FileNameList=[IndexTable(:,1)];
FileNameList=vertcat(FileNameList{:});
NumClass=14; %number of classes
JointSet=[];
% dictionary size
NumAtom=400;%---------------------------------------------------
% chunk size
nChunk=11; % for "cornell"
nChunk=13;
nFrame=5;
%sampling step
Step=1;
StepFrame=1
DataName='Cornell';         %------------------------------------------------tunable according to database
% DataName='MSRAction3D';
if strcmp(DataName,'Cornell')
nSubject=4; % for cornell database
else
nSubject=1; % for MSR Action 3D
end
LeftHand=3;
RandomAction=6;StillAction=9;
nTrainSubject=5; % for MSR Action 3D, first nTrainSubject as trainning
for j=1:nSubject
    if strcmp(DataName,'Cornell')
        TrainInd=(SubjectList~=j&ClassList~=RandomAction&ClassList~=StillAction);% remove random and still in training set   % choose training data without subject j
    else
        TrainInd=SubjectList<=nTrainSubject;
    end
    %     TrainInd=(SubjectList~=j);% remove random and still in training set   % choose training data without subject j
    BasisMatrix={};WhitenMatrix={};BasisPatch={};
    TrainClass=unique(ClassList(TrainInd));
    for k=1:length(TrainClass)
%     for k=[1:8 8:14] % remove still
        ClassInd=TrainClass(k);% index for kth class
%         ClassInd=14;% talking on the phone
        Action=ActionName{ClassInd};
        display(['training class: ' Action]);
        ClassInd=(ClassList==ClassInd);
        TrainFileName=FileNameList(ClassInd&TrainInd) % training files without jth subject
        TrainSubject=SubjectList(ClassInd&TrainInd)  % subjects for training
        data=[];Patchdata=[];
        for n=1:length(TrainFileName)
            load(fullfile(JointDatabase,[TrainFileName{n} JointFileName]));
%             if TrainSubject(n)~=LeftHand % left handed person using mirror data
                display(['get sample from non-mirror subject :' num2str(TrainSubject(n))]);  
                [data Joint]=GetSample(data,Joint,nChunk,Step,true); %sample video with chunk
%                 Patchdata=GetSample(Patchdata,PatchSample,nFrame,StepFrame,false); %sample video with chunk
%             else
%                 display(['using mirror data from subject :' num2str(TrainSubject(n))]);
%                 data=GetSample(data,MirrorJoint,nChunk,Step,true);
                
%                 Patchdata=GetSample(Patchdata,MirrorSample,nFrame,StepFrame,false); %sample video with chunk
%             end
%                  PlotVolume(Joint);
                 PlotSample(Joint,data);
        end
%         [BasisW WhitenV]=ICATrain(data);
        [BasisW WhitenV]=DictionaryLearningBasedOnSC(NumAtom,data); %learning dictionary based on SC
        title(['dictionary of ' Action]);
%         [icasig A W] = fastica(Patchdata, 'numOfIC',200);
        BasisMatrix(k)={BasisW};WhitenMatrix(k)={WhitenV}; 
%         BasisPatch(k)={W};
        Basis=(BasisW*WhitenV)';
        figure(1);plot(Basis(:,1:10:100));title(['Basis for ' ActionName{TrainClass(k)}]);
    end
    
    save([ResultBasisDir 'Without' num2str(j) 'Subject_BasisWhitenTrainIndex.mat'],'BasisMatrix','WhitenMatrix','TrainInd');
%     save([ResultBasisDir 'Without' num2str(j) 'Subject_BasisWhitenTrainIndex.mat'],'BasisMatrix','WhitenMatrix','TrainInd','BasisPatch');
% ind=find(strcmp(FileList,'activityLabel.txt')|strcmp(FileList,'README.txt'));
% FileList(ind)=[];

%%
% for i=1:length(FileList)
%     FileName=FileList{i};
%     Name=[JointDatabase FileName];
%     matrix Joint
%     load(Name); 
%     remove mean value
%     NormalizedData=Joint-repmat(mean(Joint),[size(Joint,1) 1 1]);
%     'slide' dense sampling
%     
%     NormalizedData=permute(NormalizedData,[3 2 1]);
%     nSize=size(NormalizedData);
%     data=im2colstep(NormalizedData,[nChunk nSize(2) nSize(3)],[Step 1 1]); 
%     [BasisW WhitenMatrix]=ICATrain(data);
%     save([ResultBasisDir FileName(1:end-4) '_Basis.mat'],'BasisW','WhitenMatrix');
%     for j=1:25:1000
%         close all;
%         visualizeSkeleton(Name,j);
%         Joint= read_joint(Name);
%         save([ResultDir FileName(1:end-4) 'JointPoint.mat'],'Joint');
%     end
% end
end
end
%%


function [data NormalizedData]=GetSample(data,Joint,nChunk,Step,bJoint)
%   Joint number -> Joint name
%      1 -> HEAD
%      2 -> NECK
%      3 -> TORSO
%      4 -> LEFT_SHOULDER
%      5 -> LEFT_ELBOW
%      6 -> RIGHT_SHOULDER
%      7 -> RIGHT_ELBOW
%      8 -> LEFT_HIP
%      9 -> LEFT_KNEE
%     10 -> RIGHT_HIP
%     11 -> RIGHT_KNEE
%     12 -> LEFT_HAND
%     13 -> RIGHT_HAND
%     14 -> LEFT_FOOT
%     15 -> RIGHT_FOOT

 % normalize by the distance between two shoulders
%  LeftShoulder=4;RightShoulder=6;
%  ShoulderLine=squeeze(Joint(LeftShoulder,:,:)-Joint(RightShoulder,:,:));
%  DistanceShoulder=sqrt(sum(ShoulderLine.^2));
 %%
 
 %remove mean value
 if bJoint
     NormalizedData=Joint-repmat(mean(Joint),[size(Joint,1) 1 1]);
 else % not necessary for image patches
     NormalizedData=Joint;
 end
 NormalizedData=ComputeDistanceAngle(NormalizedData);
%  NormalizedData=NormalizedData./(0.1*reshape(repmat(DistanceShoulder,size(Joint,1)*size(Joint,2),1),size(Joint,1),size(Joint,2),size(Joint,3)));
 
             % 'slide' dense sampling
 NormalizedData=permute(NormalizedData,[3 2 1]); % 1st dimension is coordiantes for the same joint to reserve continuity
 nSize=size(NormalizedData);
 if bJoint
     data=[data im2colstep(NormalizedData,[nChunk nSize(2) nSize(3)],[Step 1 1])];
 else
     data=[data im2colstep(double(NormalizedData),[nChunk nSize(2) 1],[Step 1 1])];
 end
end

function [BasisW WhitenMatrix]=ICATrain(data)
startup

%% Clear
% clear ; close all ; clc ;

%% Load Data
%  You can obtain patches.mat from 
%  http://cs.stanford.edu/~jngiam/data/patches.mat

% fprintf('Loading Data\n');
% 
% %  Loads a variable data (size 256x50000)
% load patches.mat

%  Reduce dataset size for faster training
% data = data(:, 1:20000);

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

%% Run the optimization with minFunc (ICA)
fprintf('\nTraining ICA (w/ Score Matching)\n\n');
% nHidden = 400;%original; 
nHidden = 400; 
nInput = size(whiteData, 1);
W = randn(nHidden, nInput); 
options.Method  = 'lbfgs';
options.maxIter = 100;	    % Maximum number of iterations of L-BFGS to run 
options.display = 'on';

tic
[optW, cost] = minFunc( @icaScoreMatching, ...
                        W(:), options, whiteData, ...
                        nHidden, nInput);
toc

%% Display Results
optW = reshape(optW, nHidden, nInput);
displayData(optW * V);
BasisW=optW;WhitenMatrix=V;
fprintf('ICA Training Completed.\n');
fprintf('Press Enter to Continue.\n\n');
% pause
end
function PlotSample(Joint,DataSample)
figure(1);
vol3d('CData',Joint,'Alpha',1*ones(size(Joint)));
% title('Joint Volume');
title('Distance Angle Volume');
% xlabel('xyz');ylabel('Joint');zlabel('frame');
xlabel('Distance Angle Component');ylabel('Joint');zlabel('Frame');
figure(2);
imagesc(DataSample);
title('Coordinate Sample');
xlabel('Sample Index');ylabel('coordiante xyz');
end
function DistAngle=ComputeDistanceAngle(NormalizedJoint)
nSize=size(NormalizedJoint);
ReferenceJoint=NormalizedJoint(1,:,:);
ReferenceJoint=repmat(ReferenceJoint,[size(NormalizedJoint,1),1,1]);
NormalizedJoint=mat2cell(NormalizedJoint,ones(1,nSize(1)),3,ones(1,nSize(3)));
ReferenceJoint=mat2cell(ReferenceJoint,ones(nSize(1),1),3,ones(nSize(3),1));
Distance=cellfun(@(x1,x2)norm(x1-x2),NormalizedJoint,ReferenceJoint,'UniformOutput',false);
Angle=cellfun(@(x1,x2)real(acos(dot(x1/norm(x1),x2/norm(x2)))),NormalizedJoint,ReferenceJoint,'UniformOutput',false);
DistAngle=cat(2,Distance,Angle);
DistAngle=cell2mat(DistAngle);
end
function PlotVolume(Joint)
Joint=permute(Joint,[2 1 3]);
v=Joint;
[y x z]=size(Joint);
[X Y Z]=meshgrid(1:x,1:y,1:z);
maxX=max(X(:));maxY=max(Y(:));maxZ=max(Z(:));

hx = slice(x,y,z,v,maxX,[],[]);
set(hx,'FaceColor','interp','EdgeColor','none')
end