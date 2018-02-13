function JointCorrelation()
%generate mesh
%Author: Jin Qi
%Date:   3/01/2014
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
ResultDir='./Result/';
ResultBasisDir='./Result/Basis/';
ResultIndex='./Result/Index/';
ResultCovDir='./Result/Cov/';
mkdir(ResultDir);
mkdir(ResultBasisDir);
mkdir(ResultCovDir);

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
% chunk size
nChunk=11;
nFrame=5;
%sampling step
Step=1;
StepFrame=1
nSubject=4;
LeftHand=3;
RandomAction=6;StillAction=9;
% for j=1:nSubject
%     TrainInd=(SubjectList~=j&ClassList~=RandomAction&ClassList~=StillAction);% remove random and still in training set   % choose training data without subject j
% %     TrainInd=(SubjectList~=j);% remove random and still in training set   % choose training data without subject j
%     BasisMatrix={};WhitenMatrix={};BasisPatch={};
%     TrainClass=unique(ClassList(TrainInd));
    ClassName=unique(ClassList);
    for k=1:length(ClassName)
%     for k=[1:8 8:14] % remove still
        ClassInd=ClassName(k);                   % index for kth class
        display(['training class: ' ActionName{ClassInd}]);
        ClassInd=(ClassList==ClassInd);
        FileName=FileNameList(ClassInd) % training files without jth subject
%         TrainSubject=SubjectList(ClassInd&TrainInd)  % subjects for training
%         CovClass=[];Patchdata=[];
        for n=1:length(FileName)
            load(fullfile(JointDatabase,[FileName{n} 'JointPoint.mat']));
%             if TrainSubject(n)~=LeftHand % left handed person using mirror data
%                 display(['get sample from non-mirror subject :' num2str(TrainSubject(n))]);  
                Cov_xyz=ComputeCovMatrix(Joint,nChunk,Step,true); %sample video with chunk
                CovClass(k).CovSample(n).CovMatrix=Cov_xyz;
%                 Patchdata=GetSample(Patchdata,PatchSample,nFrame,StepFrame,false); %sample video with chunk
%             else
%                 display(['using mirror data from subject :' num2str(TrainSubject(n))]);
%                 data=GetSample(data,MirrorJoint,nChunk,Step,true);
                
%                 Patchdata=GetSample(Patchdata,MirrorSample,nFrame,StepFrame,false); %sample video with chunk
%             end
        end
%         [BasisW WhitenV]=ICATrain(data);
% %         [icasig A W] = fastica(Patchdata, 'numOfIC',200);
%         BasisMatrix(k)={BasisW};WhitenMatrix(k)={WhitenV}; 
% %         BasisPatch(k)={W};
%         Basis=(BasisW*WhitenV)';
%         figure(1);plot(Basis(:,1:10:100));title(['Basis for ' ActionName{TrainClass(k)}]);
    end
    save(fullfile(ResultCovDir,'JointCovMatrix.mat'),'CovClass','ClassName','ActionName');
    PlotCovMatrix(fullfile(ResultCovDir,'JointCovMatrix.mat'));
%     save([ResultBasisDir 'Without' num2str(j) 'Subject_BasisWhitenTrainIndex.mat'],'BasisMatrix','WhitenMatrix','TrainInd');
%     save([ResultBasisDir 'Without' num2str(j) 'Subject_BasisWhitenTrainIndex.mat'],'BasisMatrix','WhitenMatrix','TrainInd','BasisPatch');
% ind=find(strcmp(FileList,'activityLabel.txt')|strcmp(FileList,'README.txt'));
% FileList(ind)=[];
% end
end
%%
function PlotCovMatrix(JointCovMatrix_FileName)
load(JointCovMatrix_FileName); % load variable "CovClass,ClassName,ActionName"
[nRow,nCol,nComp]=size(CovClass(1).CovSample(1).CovMatrix);
NumClass=length(CovClass);
nRowSubImage=ceil(sqrt(NumClass));% number of rows of subimages
for k=1:length(CovClass)
    ClassInd=ClassName(k);                   % index for kth class
    display(['training class: ' ActionName{ClassInd}]);
    KthClass=CovClass(k).CovSample;
    NumVideo=length(KthClass);
    nRowCol=ceil(sqrt(NumVideo)); % number of row and column for big matrix containing sub images
    CombinedImage=cell(nRowCol,nRowCol);
    BorderWidth=2;
    CombinedImage(:)={zeros(nRow+BorderWidth,nCol+BorderWidth)}; % 2 for border width
%     CombinedImage=zeros(nRowCol*[nRow nCol]); % big image
    for m=1:length(KthClass)
        CovMatrix=abs(sum(KthClass(m).CovMatrix,3));
        CovMatrix=padarray(CovMatrix,[BorderWidth/2 BorderWidth/2],100,'both');
        CombinedImage{ind2sub([nRowCol nRowCol],m)}=CovMatrix;
    end
    CombinedImage=cell2mat(CombinedImage);
    [R C]=ind2sub([nRowSubImage nRowSubImage],k);
    figure(1);subplot(nRowSubImage,nRowSubImage,k);
%     imshow(log(CombinedImage),[]);
    imshow((CombinedImage),[]);
%     subimage(CombinedImage,colormap(hsv));
%     imagesc(CombinedImage);
    axis square;
    title(['Action: ' ActionName{ClassInd}]);
end
end

function data=ComputeCovMatrix(Joint,nChunk,Step,bJoint)
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
%  NormalizedData=NormalizedData./(0.1*reshape(repmat(DistanceShoulder,size(Joint,1)*size(Joint,2),1),size(Joint,1),size(Joint,2),size(Joint,3)));
 
             % 'slide' dense sampling
 NormalizedData=permute(NormalizedData,[3 2 1]); % 1st dimension is coordiantes for the same joint to reserve continuity
CovArray=[];
 for i=1:size(NormalizedData,2)
     CoorComp=squeeze(NormalizedData(:,i,:)); % x,y,z coordinate component
     CovComp=cov(CoorComp);
     CovArray(:,:,i)=CovComp;
 end
 data=CovArray;
%  nSize=size(NormalizedData);
%  if bJoint
%      data=[data im2colstep(NormalizedData,[nChunk nSize(2) nSize(3)],[Step 1 1])];
%  else
%      data=[data im2colstep(double(NormalizedData),[nChunk nSize(2) 1],[Step 1 1])];
%  end
end

