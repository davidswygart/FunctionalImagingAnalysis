responseShape
% Best seems to be using dF, and computing dF for each epoch
set(groot,'defaultLineLineWidth',2.0)
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

ne = floor(size(smallGreen,4)/2)*2; % round down so that big and small have same number of epochs
%% Use best spot for both large and small SI = 0
smallG1 = smallGreen(:,:,:, 1:2:ne);
smallT1 = smallTime(:,:,:, 1:2:ne);
smallG2 = smallGreen(:,:,:, 2:2:ne);
smallT2 = smallTime(:,:,:, 2:2:ne);

%% Use 1200 preTimes as StimTimes for large spot SI = 1
bigG = bigGreen(:,:,:, 2:2:ne);
bigT = bigTime(:,:,:, 2:2:ne);

bigG(bigT > 0) = nan;
bigT(bigT > 0) = nan;
bigT = bigT * -1;

smallG3 = smallG2;
smallT3 = smallT2;

smallG3(smallT3 > 0) = nan;
smallT3(smallT3 > 0) = nan;

smallG3 = cat(3, smallG3,bigG);
smallT3 = cat(3, smallT3,bigT);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 0 suppression, Big/Small, ratio %%%%%%%%%%%%%%%%%
figure(1)
clf
hold on

% calculate suppression for each epoch and then avg them together: r = 0
[dF_SI, dFoF_SI, df_perm, dfof_perm] = calcRatio_everyEpoch_shuffled(smallG1,smallG2,smallT1,smallT2);
[f,x] = ecdf(dF_SI(meetsThresh));
f(f>.5) =  1-f(f>.5);
plot(x,f)

[f,x] = ecdf(dFoF_SI(meetsThresh));
f(f>.5) =  1-f(f>.5);
plot(x,f)

[f,x] = ecdf(df_perm(meetsThresh));
f(f>.5) =  1-f(f>.5);
plot(x,f)

[f,x] = ecdf(dfof_perm(meetsThresh));
f(f>.5) =  1-f(f>.5);
plot(x,f)

% finish plotting
xlim([-.5 .5]);
ylim([.05 .5]);
plot([0,0], [.05,.5], 'k--')
legend('df','dFoF','df Perm','dFoF perm','Location', 'south')
set(gca, 'YScale', 'log')
saveas(gcf,[image,'_Ratio0_median.png'])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot full suppression, Big/Small, ratio %%%%%%%%%%%%%%%%%
figure(2)
clf
hold on

% calculate suppression for each epoch and then avg them together: SI 1
[dF_SI, dFoF_SI, df_perm, dfof_perm] = calcRatio_everyEpoch_shuffled(smallG1,smallG3,smallT1,smallT3);
[f,x] = ecdf(dF_SI(meetsThresh));
f(f>.5) =  1-f(f>.5);
plot(x,f)

[f,x] = ecdf(dFoF_SI(meetsThresh));
f(f>.5) =  1-f(f>.5);
plot(x,f)

[f,x] = ecdf(df_perm(meetsThresh));
f(f>.5) =  1-f(f>.5);
plot(x,f)

[f,x] = ecdf(dfof_perm(meetsThresh));
f(f>.5) =  1-f(f>.5);
plot(x,f)

% finish plotting
xlim([.5 1.5]);
ylim([.05 .5]);
plot([1,1], [.05,.5], 'k--')
legend('df','dFoF','df Perm','dFoF perm','Location', 'south')
set(gca, 'YScale', 'log')
saveas(gcf,[image,'_Ratio1_median.png'])

