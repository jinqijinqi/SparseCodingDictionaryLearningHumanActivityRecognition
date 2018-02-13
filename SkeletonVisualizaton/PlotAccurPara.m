function PlotAccurPara()
%Author: Jin Qi
%Date:   8/22/2014
%Email:  jqichina@hotmail.com
%copyright2014@gru
%%
clc;
close all;
ResultPaperFigure='./Result/PaperFigure';
addpath(genpath(ResultPaperFigure));
Database='CAD60';
%% load data
FileName='AverageAccuracyPrecisionRecallParameter.txt';
fid=fopen(fullfile(ResultPaperFigure,FileName),'r');
nFramePara=15;% number of nChunk parameter value 
nWordPara=21; %number of nHidden parameter value
nCoefPara=21; %number of nCoeficient parameter value
DataCell=textscan(fid,'%f %f %f %d %d %d');
Accuracy=DataCell{1};
%%

nChunkPara=DataCell{4};
nChunkPara=nChunkPara(1:nFramePara);
hf=figure(1);set(gca,'fontsize',20);
plot(nChunkPara,Accuracy(1:nFramePara));
xlabel('Segment Length');ylabel('Average Accuracy');
title('Average Accuracy VS Segment Length');
print(hf,'-depsc2',[ResultPaperFigure '/' Database 'AverAccSegLen.eps']);

WordPara=DataCell{5};
WordPara=WordPara((nFramePara+1):(nFramePara+nWordPara));
hf=figure(1);set(gca,'fontsize',20);
plot(WordPara,Accuracy((nFramePara+1):(nFramePara+nWordPara)));
xlabel('#Number of Words');ylabel('Average Accuracy');
title('Average Accuracy VS Number of Words');
print(hf,'-depsc2',[ResultPaperFigure '/' Database 'AverAccNumWord.eps']);

CoefPara=DataCell{6};
CoefPara=CoefPara((nFramePara+nWordPara+1):(nFramePara+nWordPara+nCoefPara));
hf=figure(1);set(gca,'fontsize',15);
plot(WordPara,Accuracy((nFramePara+nWordPara+1):(nFramePara+nWordPara+nCoefPara)));
xlabel('#Number of Largest Coefficients');ylabel('Average Accuracy');
title('Average Accuracy VS Number of Largest Coefficients');
print(hf,'-depsc2',[ResultPaperFigure '/' Database 'AverAccNumLargCoef.eps']);

fclose(fid);
