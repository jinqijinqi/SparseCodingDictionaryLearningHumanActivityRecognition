function TimePerformance()
%%compute time performance
%Author: Jin Qi
%Date:   01/04/2015
%Email:  jqichina@hotmail.com
%copyright2014@CNMC
%% 
CurDir=pwd;
UtilityPath=[CurDir '/utility']; % path of Utility directory 
addpath(genpath(UtilityPath));
FullActionName={'high arm wave','horizontal arm wave','hammer','hand catch','forward punch', 'high throw', 'draw x', 'draw tick', 'draw circle', 'hand clap', 'two hand wave', 'side-boxing','bend','forward kick', 'side kick', 'jogging','tennis swing', 'tennis serve', 'golf swing', 'pickup & throw'};

AS1=[2 3 5 6 10 13 18 20];AS2=[1 4 7 8 9 11 14 12];AS3=[6 14 15 16 17 18 19 20];
ASAll=[1:20];
AS={AS1,AS2,AS3,ASAll};

TestTimePerVideo=0;
TrainTime=0;
numTrain=0;numTest=0;numClasses=0;
AveAcc=0;AvePre=0;AveRec=0;
JointDir='.\Result\Joint\MSRAction3DSkeletonReal3D\';
PaperFigureDir='./Result/PaperFigure/';

%% 

for i=1:numel(AS)
ts=tic
JointList=dir([JointDir '*.mat']);
FileList={JointList.name};
FileList=cellfun(@(x)[JointDir x],FileList,'UniformOutput',false);
cellfun(@(x)delete(x),FileList,'UniformOutput',false);

% Main_Skeleton_DeleteFalseJointFile_MSRAction3D(AS{i});
ts=tic
Main_SparseCoding_Dictionary();
tDict=toc(ts)

ts=tic
Main_HistFeature_DirectProjection_MSRAction3D();
tFeat=toc(ts)

ts=tic
conf=Main_TranTest_MSRAction3D();
tTrainTest=toc(ts)

close all;
FS=11;
PlotConfusion(conf.confus,FullActionName(AS{i}),FullActionName(AS{i}),FS);
title(sprintf('(%.2f %% accuracy)', ...
              100 * sum(diag(conf.confus)/conf.numTest) ),'fontsize',FS) ;

printpdf(gcf,[conf.ResultPaperFigure num2str(i) 'Confusion']) ;

TestTimePerVideo=TestTimePerVideo+conf.TestTimePerVideo+tFeat/conf.numTestTrain;
TrainTime=TrainTime+tDict+tFeat*conf.numTrain/conf.numTestTrain+tTrainTest-conf.TestTime;
numTrain=numTrain+conf.numTrain;numTest=numTest+conf.numTest;numClasses=numClasses+conf.numClasses;
AveAcc=AveAcc+conf.AverageAccuracy;
AvePre=AvePre+conf.AveragePrecision;
AveRec=AveRec+conf.AverageRecall;
end
%% 

n=numel(AS);
AverageTestTimePerVideo=TestTimePerVideo/n;
AverageTrainTime=TrainTime/n;
AveragenumTrain=numTrain/n;AveragenumTest=numTest/3;AveragenumClasses=numClasses/n;

AveAcc=AveAcc/n;
AvePre=AvePre/n;
AveRec=AveRec/n;
h=fopen([conf.ResultPaperFigure 'TimePerformance.txt'],'w+');
fprintf(h,'AverageTestTimePerVideo: %f\n AverageTrainTime: %f\n AveragenumTrain: %f\n AveragenumTest: %f\n AveragenumClasses: %f\n',AverageTestTimePerVideo,AverageTrainTime,AveragenumTrain,AveragenumTest,AveragenumClasses);
fclose(h);
h=fopen([conf.ResultPaperFigure 'AveAccRecPre.txt'],'w+');
fprintf(h,'Average Accuracy: %f\n average recall: %f \n average precision: %f\n',AveAcc,AveRec,AvePre);
fclose(h);