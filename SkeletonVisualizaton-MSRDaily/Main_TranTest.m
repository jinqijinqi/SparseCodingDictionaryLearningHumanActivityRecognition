function Main_TranTest()
%generate mesh
%Author: Jin Qi
%Date:   1/30/2014
%Email:  jqichina@hotmail.com
%copyright2014@gru
%%
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
ResultHistFeatureDir='./Result/HistFeature/';
ResultIndex='./Result/Index/';

load(fullfile(ResultIndex,'IndexTable.mat'));
SubjectList=cell2mat([IndexTable(:,2)]);
ClassList=cell2mat([IndexTable(:,3)]);
FileNameList=[IndexTable(:,1)];
FileNameList=vertcat(FileNameList{:});

HistFeatList=dir([ResultHistFeatureDir '*.mat']); 
HistFeatList={HistFeatList.name};
% TrainIndex=TrainInd;
%%

conf.calDir = 'data/caltech-101' ;
conf.dataDir = 'data/' ;
mkdir(conf.dataDir);
conf.autoDownloadData = false ;
conf.numTrain = 15 ;
conf.numTest =17 ;
conf.numClasses = 102 ;
conf.numWords = 600 ;
conf.numSpatialX = [2 4] ;
conf.numSpatialY = [2 4] ;
conf.quantizer = 'kdtree' ;
conf.svm.C = 10. ;

   conf.svm.solver = 'sdca' ;
%    conf.svm.solver = 'sgd' ;
  conf.svm.solver = 'liblinear' ; % 85.85 for MSRAction3D

conf.svm.biasMultiplier = 1 ;% original
 conf.svm.biasMultiplier = 0.5 ;

conf.phowOpts = {'Step', 3} ;
conf.clobber = false ;
conf.tinyProblem = true ;
conf.prefix = 'baseline' ;
conf.randSeed = 1 ;

if conf.tinyProblem
  conf.prefix = 'tiny' ;
  conf.numClasses = 5 ;
  conf.numSpatialX = 2 ;
  conf.numSpatialY = 2 ;
  conf.numWords = 300 ;
  conf.phowOpts = {'Verbose', 2, 'Sizes', 7, 'Step', 5} ;
end

conf.vocabPath = fullfile(conf.dataDir, [conf.prefix '-vocab.mat']) ;
conf.histPath = fullfile(conf.dataDir, [conf.prefix '-hists.mat']) ;
conf.modelPath = fullfile(conf.dataDir, [conf.prefix '-model.mat']) ;
conf.resultPath = fullfile(conf.dataDir, [conf.prefix '-result']) ;

randn('state',conf.randSeed) ;
rand('state',conf.randSeed) ;
vl_twister('state',conf.randSeed) ;

% --------------------------------------------------------------------
%                                            Download Caltech-101 data
% --------------------------------------------------------------------

% if ~exist(conf.calDir, 'dir') || ...
%    (~exist(fullfile(conf.calDir, 'airplanes'),'dir') && ...
%     ~exist(fullfile(conf.calDir, '101_ObjectCategories', 'airplanes')))
%   if ~conf.autoDownloadData
%     error(...
%       ['Caltech-101 data not found. ' ...
%        'Set conf.autoDownloadData=true to download the required data.']) ;
%   end
%   vl_xmkdir(conf.calDir) ;
%   calUrl = ['http://www.vision.caltech.edu/Image_Datasets/' ...
%     'Caltech101/101_ObjectCategories.tar.gz'] ;
%   fprintf('Downloading Caltech-101 data to ''%s''. This will take a while.', conf.calDir) ;
%   untar(calUrl, conf.calDir) ;
% end
% 
% if ~exist(fullfile(conf.calDir, 'airplanes'),'dir')
%   conf.calDir = fullfile(conf.calDir, '101_ObjectCategories') ;
% end

% --------------------------------------------------------------------
%                                                           Setup data
% --------------------------------------------------------------------
% classes = dir(conf.calDir) ;
% classes = classes([classes.isdir]) ;
% classes = {classes(3:conf.numClasses+2).name} ;
% 
% images = {} ;
% imageClass = {} ;
% for ci = 1:length(classes)
%   ims = dir(fullfile(conf.calDir, classes{ci}, '*.jpg'))' ;
%   ims = vl_colsubset(ims, conf.numTrain + conf.numTest) ;
%   ims = cellfun(@(x)fullfile(classes{ci},x),{ims.name},'UniformOutput',false) ;
%   images = {images{:}, ims{:}} ;
%   imageClass{end+1} = ci * ones(1,length(ims)) ;
% end
% selTrain = find(mod(0:length(images)-1, conf.numTrain+conf.numTest) < conf.numTrain) ;
% selTest = setdiff(1:length(images), selTrain) ;
% imageClass = cat(2, imageClass{:}) ;

% model.classes = classes ;
model.classes=ClassList;
model.phowOpts = conf.phowOpts ;
model.numSpatialX = conf.numSpatialX ;
model.numSpatialY = conf.numSpatialY ;
model.quantizer = conf.quantizer ;
model.vocab = [] ;
model.w = [] ;
model.b = [] ;
model.classify = @classify ;

% --------------------------------------------------------------------
%                                                     Train vocabulary
% --------------------------------------------------------------------

% if ~exist(conf.vocabPath) || conf.clobber
% 
%   % Get some PHOW descriptors to train the dictionary
%   selTrainFeats = vl_colsubset(selTrain, 30) ;
%   descrs = {} ;
%   %for ii = 1:length(selTrainFeats)
%   parfor ii = 1:length(selTrainFeats)
%     im = imread(fullfile(conf.calDir, images{selTrainFeats(ii)})) ;
%     im = standarizeImage(im) ;
%     [drop, descrs{ii}] = vl_phow(im, model.phowOpts{:}) ;
%   end
% 
%   descrs = vl_colsubset(cat(2, descrs{:}), 10e4) ;
%   descrs = single(descrs) ;
% 
%   % Quantize the descriptors to get the visual words
%   vocab = vl_kmeans(descrs, conf.numWords, 'verbose', 'algorithm', 'elkan', 'MaxNumIterations', 50) ;
%   save(conf.vocabPath, 'vocab') ;
% else
%   load(conf.vocabPath) ;
% end
% 
% model.vocab = vocab ;
% 
% if strcmp(model.quantizer, 'kdtree')
%   model.kdtree = vl_kdtreebuild(vocab) ;
% end
% 
% % --------------------------------------------------------------------
% %                                           Compute spatial histograms
% % --------------------------------------------------------------------
% 
% if ~exist(conf.histPath) || conf.clobber
%   hists = {} ;
%   parfor ii = 1:length(images)
%   % for ii = 1:length(images)
%     fprintf('Processing %s (%.2f %%)\n', images{ii}, 100 * ii / length(images)) ;
%     im = imread(fullfile(conf.calDir, images{ii})) ;
%     hists{ii} = getImageDescriptor(model, im);
%   end
% 
%   hists = cat(2, hists{:}) ;
%   save(conf.histPath, 'hists') ;
% else
%   load(conf.histPath) ;
% end
%% parameter setting
JointPoint.Gamma=0.7;
JointPoint.MirrorGamma=0.3;
JointPoint.C=10; %best value 85% for Joint point set "cornell"
JointPoint.lambda=0.02

JointPoint.Gamma=0.00001;
JointPoint.MirrorGamma=0.3;
JointPoint.C=10; %best value 89.70% for Joint point set "cornell"
JointPoint.lambda=0.00005

% JointPoint.Gamma=0.035;
% JointPoint.MirrorGamma=0.3;
% JointPoint.C=10; %best value accuracy 88.889/350 feature % for Joint point set "MSRAction3D"
% JointPoint.lambda=0.001

% JointPoint.Gamma=0.01;
% JointPoint.MirrorGamma=0.3;
% JointPoint.C=10; %best value accuracy 89.54/350  % for Joint point set "MSRAction3D"
% JointPoint.lambda=0.002

% JointAngle.Gamma=0.1;
% JointAngle.MirrorGamma=0.2;
% JointAngle.C=100000000;
% JointAngle.lambda=0.025; % best value 64%, first 100 feature, for Joint angle

% JointAngle.Gamma=0.75;
% JointAngle.MirrorGamma=0.3;
% JointAngle.C=100000000;
% JointAngle.lambda=0.025; %best value 61%, first 50 feature for Joint angle

% JointAngle.Gamma=0.5;
% JointAngle.MirrorGamma=0.6;
% JointAngle.C=100000000;
% JointAngle.lambda=0.025; %best value 66%, first 200 feature for Joint angle

Para=JointPoint;
% Para=JointAngle;
Gamma=Para.Gamma;
MirrorGamma=Para.MirrorGamma;
conf.svm.C=Para.C;
lambda=Para.lambda;

% --------------------------------------------------------------------
%                                                  Compute feature map
% --------------------------------------------------------------------
AverageAccuracy=0;AveragePrecision=0;AverageRecall=0;
for i=1:numel(HistFeatList)
HistFeat=load([ResultHistFeatureDir HistFeatList{i}]);
%%
%add random and still to training set
TrainSubject=SubjectList(HistFeat.TrainInd);
TrainSubject=unique(TrainSubject);
TrainInd=ismember(SubjectList,TrainSubject);
HistFeat.TrainInd=TrainInd;
%%

selTrain=HistFeat.TrainInd;
selTest=find(~selTrain);
selTrain=find(selTrain);
classes=unique(ClassList(HistFeat.TrainInd)); % trained classes
classes=num2cell(classes);

imageClass=ClassList';
hists=HistFeat.HistFeature;
% hists=HistFeat.HistFeaturePatch;
% hists=[hists;HistFeat.HistFeaturePatch];
psix = vl_homkermap(hists, 1, 'kchi2', 'gamma', .5) ;%original
psix = vl_homkermap(hists, 1, 'kchi2', 'gamma', Gamma) ;
%  psix=hists;

% MirrorHist=HistFeat.MirrorHistFeature;
% % MirrorHist=HistFeat.MirrorHistFeaturePatch;
% % MirrorHist=[MirrorHist;HistFeat.MirrorHistFeaturePatch];
% Mirrorpsix = vl_homkermap(MirrorHist, 1, 'kchi2', 'gamma', .5) ;%original
% Mirrorpsix = vl_homkermap(MirrorHist, 1, 'kchi2', 'gamma', MirrorGamma) ;
% Mirrorpsix=MirrorHist;

HandLeft=3;   % subject 3 is hand lefted person
HandLeftInd=(SubjectList==HandLeft);
% psix(:,HandLeftInd)=Mirrorpsix(:,HandLeftInd);
%        psix=Mirrorpsix;
% --------------------------------------------------------------------
%                                                            Train SVM
% --------------------------------------------------------------------

% if ~exist(conf.modelPath) || conf.clobber
  switch conf.svm.solver
    case {'sgd', 'sdca'}
      lambda = 1 / (conf.svm.C *  length(selTrain)) ;
%      lambda=0.001;
      w = [] ;
      for ci = 1:length(classes)
        perm = randperm(length(selTrain)) ;
        fprintf('Training model for class %s\n', classes{ci}) ;
        y = 2 * (imageClass(selTrain) == classes{ci}) - 1 ;
%         [w(:,ci) b(ci) info] = vl_svmtrain([psix(:, selTrain(perm))], [y(perm)], lambda, ...
%         [w(:,ci) b(ci) info] = vl_svmtrain([psix(:, selTrain(perm)) Mirrorpsix(:, selTrain(perm))], [y(perm) y(perm)], lambda, ...
          [w(:,ci) b(ci) info] = vl_svmtrain([psix(:, selTrain(perm))], [y(perm)], lambda, ...
          'Solver', conf.svm.solver, ...
          'MaxNumIterations', 10/lambda, ...
          'BiasMultiplier', conf.svm.biasMultiplier, ...
          'Epsilon', 1e-3);
      end

    case 'liblinear'
      for ci = 1:length(classes)
        perm = randperm(length(selTrain)) ;
        fprintf('Training model for class %s\n', classes{ci}) ;
        y = 2 * (imageClass(selTrain) == classes{ci}) - 1 ;
         [w(:,ci) b(ci) info] = vl_svmtrain([psix(:, selTrain(perm))], [y(perm)], lambda);
%       svm = vl_svmtrain(imageClass(selTrain)', ...
%                   sparse(double(psix(:,selTrain))),  ...
%                   sprintf(' -s 3 -B %f -c %f', ...
%                           conf.svm.biasMultiplier, conf.svm.C), ...
%                   'col') ;
%       w(:,ci) = svm.w(:,1:end-1)' ;
%       b(ci) =  svm.w(:,end)' ;
      end
  end

  model.b = conf.svm.biasMultiplier * b ;
  model.w = w ;
  
  save(conf.modelPath, 'model') ;
% else
%   load(conf.modelPath) ;
% end

% --------------------------------------------------------------------
%                                                Test SVM and evaluate
% --------------------------------------------------------------------

% Estimate the class of the test images
% scores = model.w' * Mirrorpsix + model.b' * ones(1,size(Mirrorpsix,2)) ;
scores = model.w' * psix + model.b' * ones(1,size(psix,2)) ;
[drop, imageEstClass] = max(scores, [], 1) ;
imageEstClass=[classes{imageEstClass}];
nTotalClass=length(unique(ClassList)); % total number of classes in all data
TempScore=zeros(nTotalClass,size(scores,2));
TempScore(cell2mat(classes),:)=scores;
scores=TempScore;
% Compute the confusion matrix
classes=unique(ClassList);
% imageClass=imageClass';
idx = sub2ind([length(classes), length(classes)], ...
              imageClass(selTest), imageEstClass(selTest)) ;
confus = zeros(length(classes)) ;
confus = vl_binsum(confus, ones(size(idx)), idx) ;

% Plots
figure(1) ; clf;
subplot(1,2,1) ;
imagesc(scores(:,[selTrain' selTest'])) ; title('Scores') ;
set(gca, 'ytick', 1:length(classes), 'yticklabel', classes) ;
subplot(1,2,2) ;
imagesc(confus) ;
conf.numTest=length(selTest);
title(sprintf('Confusion matrix (%.2f %% accuracy)', ...
              100 * sum(diag(confus)/conf.numTest) )) ;
print('-depsc2', [conf.resultPath num2str(i) '.ps']) ;
save([conf.resultPath num2str(i) '.mat'], 'confus', 'conf') ;
AverageAccuracy=AverageAccuracy+100.0 * sum(diag(confus)/conf.numTest);
AveragePrecision=AveragePrecision+100.0*mean((diag(confus))'./(sum(confus)+eps));%column sum
AverageRecall=AverageRecall+100.0*mean((diag(confus))'./(sum(confus')+eps));      %row sum
end
nTestTime=numel(HistFeatList);
xlabel(['Average Accuracy: ' num2str(AverageAccuracy/nTestTime) ' Average Precision: ' num2str(AveragePrecision/nTestTime)...
        ' Average Recall: ' num2str(AverageRecall/nTestTime)]);

end
% -------------------------------------------------------------------------
function im = standarizeImage(im)
% -------------------------------------------------------------------------

im = im2single(im) ;
if size(im,1) > 480, im = imresize(im, [480 NaN]) ; end
end
% -------------------------------------------------------------------------
function hist = getImageDescriptor(model, im)
% -------------------------------------------------------------------------

im = standarizeImage(im) ;
width = size(im,2) ;
height = size(im,1) ;
numWords = size(model.vocab, 2) ;

% get PHOW features
[frames, descrs] = vl_phow(im, model.phowOpts{:}) ;

% quantize local descriptors into visual words
switch model.quantizer
  case 'vq'
    [drop, binsa] = min(vl_alldist(model.vocab, single(descrs)), [], 1) ;
  case 'kdtree'
    binsa = double(vl_kdtreequery(model.kdtree, model.vocab, ...
                                  single(descrs), ...
                                  'MaxComparisons', 50)) ;
end

for i = 1:length(model.numSpatialX)
  binsx = vl_binsearch(linspace(1,width,model.numSpatialX(i)+1), frames(1,:)) ;
  binsy = vl_binsearch(linspace(1,height,model.numSpatialY(i)+1), frames(2,:)) ;

  % combined quantization
  bins = sub2ind([model.numSpatialY(i), model.numSpatialX(i), numWords], ...
                 binsy,binsx,binsa) ;
  hist = zeros(model.numSpatialY(i) * model.numSpatialX(i) * numWords, 1) ;
  hist = vl_binsum(hist, ones(size(bins)), bins) ;
  hists{i} = single(hist / sum(hist)) ;
end
hist = cat(1,hists{:}) ;
hist = hist / sum(hist) ;
end
% -------------------------------------------------------------------------
function [className, score] = classify(model, im)
% -------------------------------------------------------------------------

hist = getImageDescriptor(model, im) ;
psix = vl_homkermap(hist, 1, 'kchi2', 'gamma', .5) ;
scores = model.w' * psix + model.b' ;
[score, best] = max(scores) ;
className = model.classes{best} ;
end

