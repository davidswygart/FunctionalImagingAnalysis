responseShape

%Best seems to be using dF, and computing dF for each epoch
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


%% calculate dF and dFoF for each epoch and then avg them together: SI 0
[si, err] = calcSiAndErrFromDF(smallG1,smallG2,smallT1,smallT2);
scatterSIErr(si(meetsThresh), err(meetsThresh))
saveas(gcf,[image,'_siErr0.png'])

dprime_null0 = calcDprime(si(meetsThresh), err(meetsThresh));
%% calculate dF and dFoF for each epoch and then avg them together: SI 1
[si, err] = calcSiAndErrFromDF(smallG1,smallG3,smallT1,smallT3);
scatterSIErr(si(meetsThresh), err(meetsThresh))
saveas(gcf,[image,'_siErr1.png'])

dprime_null1 = calcDprime(si(meetsThresh), err(meetsThresh));
%% calculate si and error for real data
[si, err] = calcSiAndErrFromDF(smallGreen,bigGreen,smallTime,bigTime);
scatterSIErr(si(meetsThresh), err(meetsThresh))
saveas(gcf,[image,'_siErrReal.png'])

dprime_real = calcDprime(si(meetsThresh), err(meetsThresh));
%% plot cdf of dPrime

clf
hold on
cdfplot(dprime_null0)
cdfplot(dprime_null1)
cdfplot(dprime_real)

legend('n0', 'n1', 'r')

xlabel('dPrime')
ylabel('prob')
saveas(gcf,[image,'_dPrimeCompare.png'])