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
%saveas(gcf,[image,'_siErr0.png'])

si_null0 = si(meetsThresh);
si_err_null0 = err(meetsThresh);
[dprime_null0, dif_null0] = calcDprime(si(meetsThresh), err(meetsThresh));
%% calculate dF and dFoF for each epoch and then avg them together: SI 1
[si, err] = calcSiAndErrFromDF(smallG1,smallG3,smallT1,smallT3);
scatterSIErr(si(meetsThresh), err(meetsThresh))
%saveas(gcf,[image,'_siErr1.png'])

si_null1 = si(meetsThresh);
si_err_null1 = err(meetsThresh);
[dprime_null1, dif_null1]= calcDprime(si(meetsThresh), err(meetsThresh));
%% calculate si and error for real data (but use the same N)
[si, err] = calcSiAndErrFromDF(smallG1,bigGreen(:,:,:,1:2:end),smallT1,bigTime(:,:,:,1:2:end));
scatterSIErr(si(meetsThresh), err(meetsThresh))
%saveas(gcf,[image,'_siErrReal.png'])

si_real = si(meetsThresh);
si_err_real = err(meetsThresh);
[dprime_real, dif_real]= calcDprime(si(meetsThresh), err(meetsThresh));

%% plot cdf of SI values
clf
hold on
cdfplot(si_null0)
cdfplot(si_null1)
cdfplot(si_real)
xlim([-1.5 2])

[~, p]  = kstest2(si_null0,si_real);
p = num2str(round(p, 3));
m = num2str(round(mean(si_null0),3));
e = num2str(round(std(si_null0),3));
l1 = ['null0:  ', m, ' +- ', e, ' (p=', p, ')'];


[~, p]  = kstest2(si_null1,si_real);
p = num2str(round(p, 3));
m = num2str(round(mean(si_null1),3));
e = num2str(round(std(si_null1),3));
l2 = ['null1:  ', m, ' +- ', e, ' (p=', p, ')'];

m = num2str(round(mean(si_real),3));
e = num2str(round(std(si_real),3));
l3 = ['real:   ', m, ' +- ', e];
legend(l1, l2, l3, 'Location', 'northwest')

xlabel('SI')
ylabel('prob')
saveas(gcf,[image,'_SiCompare.png'])

%% plot cdf of SI differences

clf
hold on
cdfplot(dif_null0)
cdfplot(dif_null1)
cdfplot(dif_real)
xlim([0 2])

[~, p]  = kstest2(dif_null0,dif_real);
p = num2str(round(p, 3));
m = num2str(round(mean(dif_null0),3));
e = num2str(round(std(dif_null0),3));
l1 = ['null0:  ', m, ' +- ', e, ' (p=', p, ')'];


[~, p]  = kstest2(dif_null1,dif_real);
p = num2str(round(p, 3));
m = num2str(round(mean(dif_null1),3));
e = num2str(round(std(dif_null1),3));
l2 = ['null1:  ', m, ' +- ', e, ' (p=', p, ')'];

m = num2str(round(mean(dif_real),3));
e = num2str(round(std(dif_real),3));
l3 = ['real:   ', m, ' +- ', e];
legend(l1, l2, l3, 'Location', 'southeast')

xlabel('Difference in SI')
ylabel('prob')
saveas(gcf,[image,'_SiDiffCompare.png'])

%% plot cdf of dPrime
clf
hold on
cdfplot(dprime_null0)
cdfplot(dprime_null1)
cdfplot(dprime_real)

[~, p]  = kstest2(dprime_null0,dprime_real);
p = num2str(round(p, 3));
m = num2str(round(mean(dprime_null0),3));
e = num2str(round(std(dprime_null0),3));
l1 = ['null0:  ', m, ' +- ', e, ' (p=', p, ')'];


[~, p]  = kstest2(dprime_null1,dprime_real);
p = num2str(round(p, 3));
m = num2str(round(mean(dprime_null1),3));
e = num2str(round(std(dprime_null1),3));
l2 = ['null1:  ', m, ' +- ', e, ' (p=', p, ')'];

m = num2str(round(mean(dprime_real),3));
e = num2str(round(std(dprime_real),3));
l3 = ['real:   ', m, ' +- ', e];
legend(l1, l2, l3, 'Location', 'southeast')

xlabel('dPrime')
ylabel('prob')
saveas(gcf,[image,'_dPrimeCompare.png'])

xlim([0 2])