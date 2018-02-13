function TimePerformance()
%%compute time performance
%Author: Jin Qi
%Date:   1/4/2015
%Email:  jqichina@hotmail.com
%copyright2015@cnmc
%% 
%1)run the following code to reproduce confusion matrices for CAD 60 dataset
%in our IEEE THMS paper.
%2) uncomment the following code to run our program from scratch
%% 


%Main_Skeleton();
     
ts=tic
%Main_SparseCoding_Dictionary();
tDict=toc(ts)

ts=tic
%Main_HistFeature_DirectProjection();
tFeat=toc(ts)

ts=tic
conf=Main_TranTest(13,400,400); 
tTrainTest=toc(ts)

% results in the relative directory "./Results/PaperFigure/"
n=conf.NumSubset;
TestTimePerVideo=conf.TestTimePerVideo+tFeat/(conf.NumSubset*conf.numTestTrain);
TrainTime=(tDict+tFeat*conf.numTrain/conf.numTestTrain+tTrainTest-conf.TestTime*n)/n;
numTrain=conf.numTrain;numTest=conf.numTest;numClasses=conf.numClasses;

save([conf.ResultPaperFigure 'TimePerformance.mat'],'TestTimePerVideo','TrainTime','numTrain','numTest','numClasses');