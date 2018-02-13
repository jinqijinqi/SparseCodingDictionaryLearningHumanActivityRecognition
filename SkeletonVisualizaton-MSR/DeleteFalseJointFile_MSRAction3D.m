function DeleteFalseJointFile_MSRAction3D()
%generate mesh
%Author: Jin Qi
%Date:   3/26/2014
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
ImageDatabase='./Database/';
ResultDir='./Result/Joint/';
mkdir(ResultDir);
ResultIndex='./Result/Index/';
mkdir(ResultIndex);
SubjectList=dir(ImageDatabase);
SubjectList=SubjectList(3:end);
SubjectNameList={SubjectList.name};
NameActivityLabel='activityLabel.txt';
IndexTable={};
nSubject=length(SubjectNameList);
nJointPerFrame=20; %number of joints per frame

TrueJointFile = textread(fullfile(TrueJointDir,'TrueJoint.txt'), '%s', 'delimiter', '\n', ...
                'whitespace', '');
TrueJointFile=cellfun(@(x)[x '_skeleton3D.txt'],TrueJointFile,'UniformOutput',false);

for j=1:length(SubjectNameList)
    close all;
    SubjectName=SubjectNameList{j};
    mkdir([ResultDir,SubjectName]);
    %     List=dir(fullfile(ImageDatabase,SubjectName,'*.txt'));
    %     FileList={List.name};
    %%
    
    % read activity label
    FileList=dir(fullfile(ImageDatabase,SubjectName,'*.txt')); 
    FileList={FileList.name};
    FileList=FileList(ismember(FileList,TrueJointFile));% remove false joint files
    FileList=strrep(FileList,'.txt','');
    [Action R]=strtok(FileList,'_');
    [Subject R]=strtok(R,'_');
    
    Action=strrep(Action,'a','');
    Subject=strrep(Subject,'s','');
    Action=cellfun(@(x)str2num(x),Action);
    Subject=cellfun(@(x)str2num(x),Subject);
    
%     ActivityData=importdata(fullfile(ImageDatabase,SubjectName,NameActivityLabel));
%     ActivityData=ActivityData(1:end-1);     % remove 'END'
%     [FileList R]= strtok(ActivityData,','); % extract file names 
%     [Action R]=strtok(R,',');               % extract action names
    [ActionName ia ic]=unique(Action);      % ic is class
    ActionName=mat2cell(ActionName,1,ones(size(ActionName,2),1));
    %%
    
    File2ClassMap{j}=containers.Map(FileList,ic);
    
    FileList=cellfun(@(x)[SubjectName '/' x],FileList,'UniformOutput', false);
    SubjectFile2ClassMap{j}=containers.Map(FileList,ic);
%     IndexTable(end+1,:)={FileList,j*ones(length(FileList),1),ic};
    IndexTable(end+1,:)={FileList',Subject',ic'};
    display('Sorted actions:\n');
    ActionName
%     ind=find(strcmp(FileList,'activityLabel.txt')|strcmp(FileList,'README.txt'));
%     FileList(ind)=[];
    for i=1:length(FileList)
        FileName=[FileList{i} '.txt'];
%         FileName=[FileList{i}];
        Name=fullfile(ImageDatabase,FileName);
%     for j=1:25:1000
%         close all;
%         visualizeSkeleton(Name,j);
        Joint=importdata(Name);
        Joint=Joint(:,1:end-1);
        Joint=mat2cell(Joint,ones(size(Joint,1),1),size(Joint,2));
        Joint=reshape(Joint,nJointPerFrame,1,[]);
        Joint=cell2mat(Joint);
%         Joint= read_joint(Name);
%         MirrorJoint=MirrorJointData(Joint);
%         PatchSample=Joint2Pixel(Joint,Name);
%         PatchSample=reshape(cell2mat(PatchSample),size(PatchSample{1},1),size(PatchSample{1},2),numel(PatchSample));
%         
%         MirrorSample=MirrorPatch(PatchSample);
%         figure(1);title(['class ' ActionName(ic(i))]);
%         PlotJointCurve(Joint);
        ClassID=ic(i);
        save(fullfile(ResultDir,[FileName(1:end-4) 'JointPoint.mat']),'Joint','ClassID');
%         save(fullfile(ResultDir,[FileName(1:end-4) 'JointPoint.mat']),'Joint','MirrorJoint','ClassID');
%         save(fullfile(ResultDir,[FileName(1:end-4) 'JointPoint.mat']),'Joint','MirrorJoint','ClassID','PatchSample','MirrorSample');
%     end
    end
end
save(fullfile(ResultIndex,'File2ClassMap.mat'),'File2ClassMap');
save(fullfile(ResultIndex,'SubjectFile2ClassMap.mat'),'SubjectFile2ClassMap','ActionName');
save(fullfile(ResultIndex,'IndexTable.mat'),'IndexTable','ActionName');
end

function SamplePatch=MirrorPatch(SamplePatch)
SwapIndex=[1 2 3 6 7 4 5 10 11 8 9 13 12 15 14];
SamplePatch=SamplePatch(SwapIndex,:,:);
end

function PatchSample=Joint2Pixel(Joint,Name)
% Here are two useful equation that can be used to convert (x,y,z) coordinate into (x,y) pixel
% number in 2D image approximately. Given the (x,y,z) location of the joint, you can find where
% the joint is located within RGBD images.
% 
% x = 156.8584456124928 + 0.0976862095248 * x - 0.0006444357104 * y + 0.0015715946682 * z
% y = 125.5357201011431 + 0.0002153447766 * x - 0.1184874093530 * y - 0.0022134485957 * z
px=[0.0976862095248 -0.0006444357104  0.0015715946682]';
py=[0.0002153447766 -0.1184874093530 -0.0022134485957]';
xc=156.8584456124928;yc=125.5357201011431;
[m n s]=size(Joint);
CellJoint=squeeze(mat2cell(Joint,m,n,ones(s,1)));
JointInImage=cellfun(@(x)x*[px py]+repmat([xc yc],size(x,1),1),CellJoint,'UniformOutput',false);

VideoDirectory=Name(1:end-4); % video directory
% FileList=new_dir(fullfile(VideoDirectory,'Depth*.png'));
FileList=new_dir(fullfile(VideoDirectory,'RGB*.png'));
FileNameList={FileList.name};
FullFileNameList=cellfun(@(x)fullfile(VideoDirectory,x),FileNameList,'UniformOutput',false);

PatchSample=cellfun(@CropIm,FullFileNameList',JointInImage,'UniformOutput',false);
end


function PatchSample=CropIm(FileName,JointFrame)
I=imread(FileName);
if size(I,3)>1
    I=rgb2gray(I);
end
[m n]=size(I);

BlockSize=21;
HalfBlock=floor(BlockSize/2);
I=padarray(I(:,:,1),size(I));
JointFrame=JointFrame+repmat([n m],size(JointFrame,1),1)-HalfBlock;
PatchSample=[];
% figure(1);imshow(I,[]);hold on; plot(JointFrame(:,1),JointFrame(:,2),'r.');
for i=1:size(JointFrame,1)
    Jx=JointFrame(i,1);Jy=JointFrame(i,2);
    rect=[Jx,Jy,BlockSize,BlockSize];
    Sample=imcrop(I,rect);
%     hold on;rectangle('position',rect);
    PatchSample=[PatchSample;Sample(:)'];
end


end

function PlotJointCurve(Joint)
Color={'y' 'm' 'c' 'r' 'g' 'b' 'k'};
MyMap=rand(size(Joint,1),3);
nJoint=size(Joint,1);
for i=1:size(Joint,1);
%     figure(1);hold on;plot3(squeeze(Joint(i,1,:)),squeeze(Joint(i,2,:)),squeeze(Joint(i,3,:)),'-','Color',MyMap(i,:));
    ColorIndex=mod(i,numel(Color));
    if ColorIndex==0
        ColorIndex=numel(Color);
    end
    figure(1);hold on;plot3(squeeze(Joint(i,1,:)),squeeze(Joint(i,2,:)),squeeze(Joint(i,3,:)),'-','color',Color{ColorIndex});
end
end

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
function MirrorJoint=MirrorJointData(Joint)
    nJoint=size(Joint,1);
    SwapIndex=[1 2 3 6 7 4 5 10 11 8 9 13 12 15 14];
%     MirrorJoint=Joint(SwapIndex,:,:);
    LeftHand=12;RightHand=13;
    Torso=3;LeftShoulder=4;RightShoulder=6;LeftHip=8,RightHip=10;% fit these 5 points to find the plane which torso reside
    ShoulderLine=squeeze(Joint(LeftShoulder,:,:)-Joint(RightShoulder,:,:));
    HipLine=squeeze(Joint(LeftHip,:,:)-Joint(RightHip,:,:));
    ShoulderLine=ShoulderLine/norm(ShoulderLine);
    HipLine=HipLine/norm(HipLine);
    Coeff=squeeze(0.5*(ShoulderLine+HipLine));
    Coeff=mat2cell(Coeff,size(Coeff,1),ones(size(Coeff,2),1));
    
%     Index=[Torso LeftShoulder RightShoulder LeftHip RightHip];
%     X=Joint(Index,1,:);Y=Joint(Index,2,:);Z=Joint(Index,3,:);
%     X=squeeze(X);Y=squeeze(Y);Z=squeeze(Z);
%     X=mat2cell(X,size(X,1),ones(size(X,2),1));Y=mat2cell(Y,size(Y,1),ones(size(Y,2),1));Z=mat2cell(Z,size(Z,1),ones(size(Z,2),1));

%     SurfFit=cellfun(@(x,y,z)fit([x y],z,'poly11'),X,Y,Z,'UniformOutput',false);
%     Coeff=cellfun(@(x)coeffvalues(x),SurfFit,'UniformOutput',false);
%     Coeff=cellfun(@FirstEntry2One,Coeff,'UniformOutput',false);
%     Coeff=cellfun(@(x)circshift(x,[0 -1]),Coeff,'UniformOutput',false); % plane of body
    
    %     X=vertcat(X{:}); Y=vertcat(Y{:}); Z=vertcat(Z{:});
    
%     SurfFit=fit([PointArray(1,index)' PointArray(2,index)'],PointArray(3,index)','poly11');
    
%     Coeff(1)=1;Coeff=circshift(Coeff,[0 -1]);
    
%     Coeff=lsqlin([x' y' z' ones(length(x),1)],zeros(length(x),1),[],[],[1 1 1 0],1);% constrained regression with a+b+c=1;
%     NormDirect=Coeff(1:3);
    NormDirect=cellfun(@(x)x(1:3),Coeff,'UniformOutput',false);
%     NormDirect=NormDirect/norm(NormDirect);%normalize
    NormDirect=cellfun(@(x)x/norm(x),NormDirect,'UniformOutput',false);%normalize
    
%     NormDirect=vertcat(NormDirect{:});
    TorsoPoint=squeeze(Joint(Torso,:,:));
    TorsoPoint=mat2cell(TorsoPoint,3,ones(size(TorsoPoint,2),1));
    
    d=cellfun(@(x,y)-dot(x',y),NormDirect,TorsoPoint,'UniformOutput',false);
    
%     NormDirect=cellfun(@(x)repmat(x,nJoint,1),NormDirect,'UniformOutput',false);
%     d=cellfun(@(x)repmat(x,nJoint,1),d,'UniformOutput',false);
    CellJoint=squeeze(mat2cell(Joint,size(Joint,1),size(Joint,2),ones(size(Joint,3),1)));
    CellJoint=CellJoint';
    
    MirrorJoint=cellfun(@(p,n,d)p-kron(2*(p*n+d'),n'),CellJoint,NormDirect,d,'UniformOutput',false);
    MirrorJoint=horzcat(MirrorJoint{:});
    MirrorJoint=reshape(MirrorJoint,size(Joint,1),size(Joint,2),size(Joint,3));
    
    MirrorJoint=MirrorJoint(SwapIndex,:,:);
%     FinalNorm(:,end+1)=NormDirect;
%     MirrorJoint=Joint;
%     MirrorJoint(:,3,:)=-1*Joint(:,3,:); % 3rd column is Z coordinate, reflect wrt x-y plane;
%     MirrorJoint=-1*MirrorJoint; %
%% draw for debug

%     figure(2);title('left handed to right handed');
%     ShoulderHipPoint=squeeze(Joint([LeftShoulder RightShoulder LeftHip RightHip],:,1));
%     hold on;plot3(ShoulderHipPoint(:,1),ShoulderHipPoint(:,2),ShoulderHipPoint(:,3),'r-');
%     FirstTorso=TorsoPoint{1};
%     NormalEndPoint=FirstTorso+100*NormDirect{1};
%     NormalLine=[FirstTorso';NormalEndPoint'];
%     hold on;plot3(NormalLine(:,1), NormalLine(:,2), NormalLine(:,3),'g-');
%     MirrorPair=[squeeze(Joint(LeftHand,:,1));squeeze(MirrorJoint(RightHand,:,1))];
%     hold on;scatter3(MirrorPair(:,1),MirrorPair(:,2),MirrorPair(:,3),'c*');
end
function x=FirstEntry2One(x)
x(1)=-1;
end
