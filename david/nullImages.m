responseShape

% Best seems to be using dF, and computing dF for each epoch
%% make threshold using only the good times
pixFrac = .1;
isGood = raw.time > minMaxTime(1) & raw.time < minMaxTime(2);
goodThresh = raw.thresh;
goodThresh(~isGood) = nan;
[meetsThresh, tProj] = thresholdFromTProjection(goodThresh, pixFrac);

%% Choose the best size
smallGreen = greenSpots{bestInd};
smallTime = timeSpots{bestInd};
bigGreen = greenSpots{end};
bigTime = timeSpots{end};

%% Use best spot for both large and small SI = 0
smallG1 = smallGreen(:,:,:, 1:2:end);
smallT1 = smallTime(:,:,:, 1:2:end);
smallG2 = smallGreen(:,:,:, 2:2:end);
smallT2 = smallTime(:,:,:, 2:2:end);

%% Use 1200 preTimes as StimTimes for large spot SI = 1
ne = size(smallG2,4);
bigG = bigGreen(:,:,:, 2:2:ne*2);
bigT = bigTime(:,:,:, 2:2:ne*2);

bigG(bigT > 0) = nan;
bigT(bigT > 0) = nan;
bigT = bigT * -1;

smallG3 = smallG2;
smallT3 = smallT2;

smallG3(smallT3 > 0) = nan;
smallT3(smallT3 > 0) = nan;

smallG3 = cat(3, smallG3,bigG);
smallT3 = cat(3, smallT3,bigT);

%% calculate using the raw values averaged accross all epochs: SI 0
[dF_SI, dFoF_SI] = calcSiRaw(smallG1,smallG2,smallT1,smallT2);
plotSIHist(dF_SI(meetsThresh), dFoF_SI(meetsThresh), 0)
saveas(gcf,[image,'_SI0_rawVals.png'])

%% calculate using the pre and stim averaged accross for each epoch: SI 0
[dF_SI, dFoF_SI] = calcSiPreStimEpochs(smallG1,smallG2,smallT1,smallT2);
plotSIHist(dF_SI(meetsThresh), dFoF_SI(meetsThresh), 0)
saveas(gcf,[image,'_SI0_preStim.png'])

%% calculate dF and dFoF for each epoch and then avg them together: SI 0
[dF_SI, dFoF_SI] = calcSi_dFEpochs(smallG1,smallG2,smallT1,smallT2);
plotSIHist(dF_SI(meetsThresh), dFoF_SI(meetsThresh), 0)
saveas(gcf,[image,'_SI0_avgDF.png'])

% %% calculate dF and dFoF for each epoch and then avg them together: SI 0
% [dF_SI, dFoF_SI] = calcSi_everyEpoch(smallG1,smallG2,smallT1,smallT2);
% plotSIHist(dF_SI(meetsThresh), dFoF_SI(meetsThresh), 0)
% saveas(gcf,[image,'_SI0_avgSI.png'])





%% calculate using the raw values averaged accross all epochs: SI 1
[dF_SI, dFoF_SI] = calcSiRaw(smallG1,smallG3,smallT1,smallT3);
plotSIHist(dF_SI(meetsThresh), dFoF_SI(meetsThresh), 1)
saveas(gcf,[image,'_SI1_rawVals.png'])

%% calculate using the pre and stim averaged accross for each epoch: SI 1
[dF_SI, dFoF_SI] = calcSiPreStimEpochs(smallG1,smallG3,smallT1,smallT3);
plotSIHist(dF_SI(meetsThresh), dFoF_SI(meetsThresh), 1)
saveas(gcf,[image,'_SI1_preStim.png'])

%% calculate dF and dFoF for each epoch and then avg them together: SI 1
[dF_SI, dFoF_SI] = calcSi_dFEpochs(smallG1,smallG3,smallT1,smallT3);
plotSIHist(dF_SI(meetsThresh), dFoF_SI(meetsThresh), 1)
saveas(gcf,[image,'_SI1_avgDF.png'])

% %% calculate dF and dFoF for each epoch and then avg them together: SI 1
% [dF_SI, dFoF_SI] = calcSi_everyEpoch(smallG1,smallG3,smallT1,smallT3);
% plotSIHist(dF_SI(meetsThresh), dFoF_SI(meetsThresh), 1)
% saveas(gcf,[image,'_SI1_avgSI.png'])


%% show N
ne