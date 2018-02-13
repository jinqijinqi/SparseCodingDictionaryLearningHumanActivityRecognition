function Main_HistFeature_SC(nChunk,FirstN)
%generate mesh
%Author: Jin Qi
%Date:   11/24/2014
%Email:  jqichina@hotmail.com
%copyright2014@CNMC
%%  
% clc;
% clear;
global PathAdded;% global variable for adding path
PathAdded=false;

CurDir=pwd;
UtilityPath=[CurDir '/utility']; % path of Utility directory 
addpath(genpath(UtilityPath));
PathAdded=true;
%%
JointDatabase='./Result/Joint/';
% JointDatabase='./Result/JointAngle/';

ResultDir='./Result/';
ResultBasisDir='./Result/Basis/';
ResultHistFeatureDir='./Result/HistFeature/';
mkdir(ResultHistFeatureDir);
ResultIndex='./Result/Index/';
% joint position
% ListJoint=dir([JointDatabase '*.mat']);
% FileListJoint={ListJoint.name};
%dictionary
% ListBasis=dir([ResultBasisDir '*.mat']);
% FileListBasis={ListBasis.name};
% ind=find(strcmp(FileList,'activityLabel.txt')|strcmp(FileList,'README.txt'));
% FileList(ind)=[];

% chunk size
if ~exist('nChunk','var')
nChunk=11; %for "cornell"
end
nChunk=13; %for "MSRACtion3D" 89.54
nFrame=5;
%sampling step
Step=1;
StepFrame=1

%%
load(fullfile(ResultIndex,'IndexTable.mat'));
SubjectList=cell2mat([IndexTable(:,2)]);
ClassList=cell2mat([IndexTable(:,3)]);
FileNameList=[IndexTable(:,1)];
FileNameList=vertcat(FileNameList{:});

BasisList=dir([ResultBasisDir '*.mat']); 
BasisList={BasisList.name};
FullData={}; % all sample data
for i=1:length(FileNameList)
%         FileNameJoint=FileListJoint{i};
%         Name=[JointDatabase FileNameJoint];
        % matrix Joint
        load([JointDatabase FileNameList{i} 'JointPoint.mat']); 
        %remove mean value
        Data=GetSample(Joint,nChunk,Step,true);
        Data=NormalizeMatrix(Data);
        
        MirrorData=GetSample(MirrorJoint,nChunk,Step,true);
        MirrorData=NormalizeMatrix(MirrorData);
%         PatchData=GetSample(PatchSample,nFrame,StepFrame,false);
%         MirrorPatchData=GetSample(MirrorSample,nFrame,StepFrame,false);
        
%         FullData(end+1)={{Data,MirrorData,PatchData,MirrorPatchData}};% nested cell array
        FullData(end+1)={{Data,MirrorData}};% nested cell array,  "cornell"
%          FullData(end+1)={{Data}};% nested cell array, "MSRAction 3D
end
for k=1:length(BasisList)
    load([ResultBasisDir BasisList{k}]); %load ICA basis and whiten matrix and training index,variables:BasisMatrix,WhitenMatrix,TrainIndex
    %normalize basis
%     BasisMatrix=cellfun(@(x)(x'*1.0/norm(x'))',BasisMatrix,'UniformOutput',false);
%     BasisMatrix=BasisMatrix';

    BasisMatrix=NormalizeBasis(BasisMatrix,WhitenMatrix);
    
    HistFeature=[];
    MirrorHistFeature=[];
    HistFeaturePatch=[];
    MirrorHistFeaturePatch=[];
    for m=1:numel(FullData)              % get histogram feature foe each video
        ActionName{ClassList(m)}
        Data=FullData{m};
        SingleVideo=Data{1};
        if exist('FirstN','var')
        HistFeat=ProjectData(SingleVideo,BasisMatrix,WhitenMatrix,FirstN);
        else
            HistFeat=ProjectData(SingleVideo,BasisMatrix,WhitenMatrix);
        end
        HistFeature=[HistFeature HistFeat];
        
        MirrorVideo=Data{2};
        if exist('FirstN','var')
        HistFeat=ProjectData(MirrorVideo,BasisMatrix,WhitenMatrix,FirstN);
        else
            HistFeat=ProjectData(MirrorVideo,BasisMatrix,WhitenMatrix);
        end
        MirrorHistFeature=[MirrorHistFeature HistFeat];
        
        %Compute image patch features
%         load([JointDatabase FileNameList{i} 'JointPoint.mat']);
%         PatchData=GetSample(PatchSample,nFrame,StepFrame,false);
%         MirrorPatchData=GetSample(MirrorSample,nFrame,StepFrame,false);
%         
%         HistFeat=ProjectData(PatchData,BasisPatch);
%         HistFeaturePatch=[HistFeaturePatch HistFeat];
%         
%         HistFeat=ProjectData(MirrorPatchData,BasisPatch);
%         MirrorHistFeaturePatch=[MirrorHistFeaturePatch HistFeat];
        
        figure(1);
%         plot(HistFeature);title(['Features without seeing subject: ' num2str(k)]);
%         plot(HistFeat);title(['Features without seeing subject: ' num2str(k)]);
        ratio=sum(HistFeat~=0)/length(HistFeat);
        plot(HistFeat);title(['Sparse histogram of "' ActionName{ClassList(m)} '" video (nonzero ratio: ' num2str(ratio) ')']);
%         figure(2);
%         plot(MirrorHistFeature);title(['Mirror Features without seeing subject: ' num2str(k)]);
    end   
    BasisName=BasisList{k};
%     save([ResultHistFeatureDir BasisName(1:end-4) '_HistFeat.mat'],...
%         'HistFeature','TrainInd'); %with train index, "MSRAction 3D"
    save([ResultHistFeatureDir BasisName(1:end-4) '_HistFeat.mat'],...
        'HistFeature','MirrorHistFeature','TrainInd'); %with train index
%     save([ResultHistFeatureDir BasisName(1:end-4) '_HistFeat.mat'],...
%         'HistFeature','MirrorHistFeature','TrainInd','HistFeaturePatch','MirrorHistFeaturePatch'); %with train index
%     CatchHist=zeros(size(BasisMatrix,1)*length(FileListBasis),1);
%     
%     
%     for j=1:length(FileListBasis)
%         FileNameBasis=FileListBasis{i};
%         load([ResultBasisDir FileNameBasis]);
%         WhitenV=WhitenMatrix;
%         WhiteData=WhitenData(data,WhitenV);
%         ProjectCoeff=[ProjectCoeff;BasisW*WhiteData];
%     end
%     [maxV ind]=max(ProjectCoeff);
%     Table=tabulate(ind);
%     CatchHist=zeros(size(BasisW,1)*length(FileListBasis),1);
%     indx=round(Table(:,1));
%     %sum(CatchHist) is 1
%     CatchHist(indx)=CatchHist(indx)+Table(:,3);
%     CatchHist=[CatchHist HistFeature];
    
%     for j=1:25:1000
%         visualizeSkeleton(Name,j);
%         Joint= read_joint(Name);
%         save([ResultDir FileName(1:end-4) 'JointPoint.mat'],'Joint');
%     end
end
end
%%
function Basis=NormalizeBasis(BasisMatrix,WhitenMatrix)
Basis=cellfun(@(x,y)x*y,BasisMatrix,WhitenMatrix,'UniformOutput',false);
% Basis=cellfun(@(x)x./norm(x),Basis,'UniformOutput',false);
Basis=cell2mat(Basis');
end

function NormalizedData=NormalizeMatrix(DataMatrix)
[m n]=size(DataMatrix);
CellMatrix=mat2cell(DataMatrix,m,ones(n,1));
NormalizedMatrix=cellfun(@(x)x./norm(x),CellMatrix,'UniformOutput',false);
NormalizedData=cell2mat(NormalizedMatrix);
end


function HistFeat=ProjectData(SingleVideo,BasisMatrix,WhitenMatrix,FirstN)
% parameters for sparse coding
noise_var= 0.01;%orginal
% noise_var=0.001;
beta= 2.2;%orginal
% beta=0.2
sigma=0.316; %orginal
% sigma=0.03
tol=.01; %orginal
Thr=0.02; %threshold for sparse coefficients
%% 
% FirstN=20;% best for joint point "Cornell" 71 for MSRAction3D
% FirstN=100;% for MSRAcion3D 82.49
% FirstN=200;% for MSRAcion3D  accuracy 86.5
% FirstN=300;% for MSRAcion3D  accuracy 87.87
% FirstN=350;% for MSRAcion3D  accuracy 89.54, accuracy 89.70 for "cornell"
if ~exist('FirstN','var')
 FirstN=400;% for MSRAcion3D  accuracy 88.89, accracy 91.17 for "cornell"
end
% FirstN=450;% 
% FirstN=200;%for joint angle
if exist('WhitenMatrix')
%     WhitenedData=cellfun(@(x)x*SingleVideo,WhitenMatrix,'UniformOutput',false);
%     ProjectCoeff=cellfun(@(x,y)x*y,BasisMatrix,WhitenedData,'UniformOutput',false);
%     ProjectCoeff=BasisMatrix*SingleVideo;   %orginal
%     ProjectCoeff=cgf_fitS(BasisMatrix',SingleVideo,noise_var,beta,sigma,tol);%slow(1995)
     BasisMatrix=NormalizeMatrix(BasisMatrix');
     ProjectCoeff=SCOMP(BasisMatrix,SingleVideo); %faster OMP
else
    ProjectCoeff=cellfun(@(x,y)x*SingleVideo,BasisMatrix,'UniformOutput',false);
end
%         ProjectCoeff=ProjectCoeff';
% %         ProjectCoeff(9)=[];% remove still action basis
% %         ProjectCoeff{6}=[]; % remove random action
%         ProjectCoeff=vertcat(ProjectCoeff{:});
          % the following code has some bug in histogram compute 
%         [maxV ind]=max(abs(ProjectCoeff)); % maximize the abosolute values
%         [maxV ind]=max((ProjectCoeff)); % maximize the abosolute values
%         [maxV ind]=sort(abs(ProjectCoeff)>thr,'descend'); % maximize the abosolute values
%         ind=ind(1:FirstN,:);
%         Table=tabulate(ind(:));
%         indx=round(Table(:,1));
%         HistFeat=zeros(size(ProjectCoeff,1),1);
%         HistFeat(indx)=Table(:,3);
        ProjectCoeff=abs(ProjectCoeff);
        HistFeat=GetHist(ProjectCoeff,Thr);
        %average 
%         HistFeat=mean(ProjectCoeff');
end

function HistFeat=GetHist(ProjectCoeff,Thr)
[r c]=find(ProjectCoeff>=Thr); % ignore small coefficients
Table=tabulate(r(:));
HistFeat=zeros(size(ProjectCoeff,1),1);
HistFeat(Table(:,1))=Table(:,3);
end

function data=GetSample(Joint,nChunk,Step,bJoint)
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
 nSize=size(NormalizedData);
 if bJoint
     data=[im2colstep(NormalizedData,[nChunk nSize(2) nSize(3)],[Step 1 1])];
 else
     data=[im2colstep(double(NormalizedData),[nChunk nSize(2) 1],[Step 1 1])];
 end
end

function WhiteData=WhitenData(data,WhitenV)
%% PCA Whitening
fprintf('\n preprocess data\n');

% Remove DC
data = bsxfun(@minus, data, mean(data, 1));

% Remove the "mean" patch
data = bsxfun(@minus, data, mean(data, 2));


%% Whiten the data
WhiteData = WhitenV * data;
end
