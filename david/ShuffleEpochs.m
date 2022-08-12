set(groot,'defaultLineLineWidth',2.0)

responseShape

%% Choose the best size
smallGreen = greenSpots{bestInd};
smallTime = timeSpots{bestInd};
%bigGreen = greenSpots{end};
%bigTime = timeSpots{end};

%%%%%%%%%%%%%%%%%% SI == 0 %%%%%%%%%%%%%%%%
bigGreen = smallGreen(:,:,:,2:2:end);
bigTime = smallTime(:,:,:,2:2:end);
smallGreen = smallGreen(:,:,:,1:2:end);
smallTime = smallTime(:,:,:,1:2:end);
%%%%%%%%%%%%%%%%%% SI == 0 %%%%%%%%%%%%%%%%

ns = size(smallGreen,4);
nb = size(bigGreen,4);
ne = min(ns,nb);

smallGreen = smallGreen(:,:,:,1:ne);
smallTime = smallTime(:,:,:,1:ne);
bigGreen = bigGreen(:,:,:,1:ne);
bigTime = bigTime(:,:,:,1:ne);


%% make threshold using SNR cutoff or intensity
figure(1)
clf
[meetsThresh, tProj] = thresholdBySNR(smallGreen,smallTime,.5);

% make threshold using intensity
% isGood = raw.time > minMaxTime(1) & raw.time < minMaxTime(2);
% goodThresh = raw.thresh;
% goodThresh(~isGood) = nan;
% [meetsThresh, tProj] = thresholdFromTProjection(goodThresh, .1);

saveas(gcf,[image,'_threshold.png'])
%%
if sum(meetsThresh(:)) > 4
%% calc suppression for each epoch (using ratio)
si = calcRatioEachEpoch(smallGreen, bigGreen, smallTime, bigTime);


%allSI = median(si,3)

%% extract the important pixels
[rs,cs] = find(meetsThresh);
pxSis = nan(size(si,3),size(rs,1));
for i = 1:size(rs,1)
    pxSis(:,i) = squeeze(si(rs(i),cs(i),:));
end

%% calc median si
realSi = median(pxSis,1);
%realSi = mean(pxSis,1);

%% shuffle and recalc
[ne, np] = size(pxSis);

if np < 100 
    numShuffles = 1E5;
else
    numShuffles = round(1E9 / np^2);
end

allShuf = nan(numShuffles,np);
for i = 1:numShuffles
    [~, pxPerms] = sort(rand(ne, np),2);
    shufSi = pxSis(pxPerms);
    allShuf(i,:) = median(shufSi,1);
    %allShuf(i,:) = mean(shufSi,1);
end
%% plot si values cdf
figure(2)
clf
[f,x] = ecdf(realSi);
f(f>.5) =  1-f(f>.5);
plot(x,f)
hold on


newX = linspace(mean(min(allShuf)),mean(max(allShuf)), 200);
newF = nan(numShuffles, length(newX));

for i = 1:numShuffles
[f,x] = ecdf(allShuf(i,:));
x(1) = x(1)-.0001;%needs to be unique
newF(i,:) = interp1(x,f,newX,'nearest','extrap');
end

avgF = mean(newF,1);
errF = std(newF,0,1);
avgF(avgF>.5) =  1-avgF(avgF>.5);
errorbar(newX,avgF,errF)



%legend('real','shuffled')
xlabel('suppression')
%xlim([-1,2])
ylim([0 .5]);
set(gca, 'YScale', 'linear')


plot([0,0],ylim, '--k')
plot([1,1],ylim, '--k')
saveas(gcf,[image,'_realSi.png'])

% %% examine kurtosis
% kurtReal = kurtosis(realSi,0);
% kurtShuf = sort(kurtosis(allShuf,0,2));
% 
% figure(3)
% clf
% histogram(kurtShuf)
% hold on
% plot([kurtReal,kurtReal],ylim)
% p = 1 - mean(kurtShuf < kurtReal);
% text(kurtReal,(max(ylim)/2), [' p = ', num2str(p)])
% xlabel('kurtosis')
% saveas(gcf,[image,'_kurtosis.png'])
% 
% 
% %% find the distance between these pixels
% dist = pdist2([rs,cs], [rs,cs]);
% isU = logical(tril(ones(size(dist)),-1));
% dist = dist(isU);
% 
% %% examine pairwise differences in suppression
% realDiff = abs(realSi - realSi');
% realDiff = realDiff(isU);
% 
% 
% figure(4)
% clf
% hold on
% cdfplot(realDiff);
% 
% fakeDiff = nan(numShuffles, length(realDiff));
% 
% fakeDiffCDF_X = linspace(mean(min(realDiff)),mean(max(realDiff)), 200);
% fakeDiffCDF_Y = nan(numShuffles, length(fakeDiffCDF_X));
% 
% for i = 1:numShuffles
%     d = abs(allShuf(i,:) - allShuf(i,:)');
%     d = d(isU);
%     fakeDiff(i,:) = d;
%     
%     [f,x] = ecdf(d);
%     x(1) = x(1)-.0001;%needs to be unique
%     fakeDiffCDF_Y(i,:) = interp1(x,f,fakeDiffCDF_X,'nearest','extrap');
% end
% 
% avgFakeDiffCDF = mean(fakeDiffCDF_Y,1);
% errFakeDiffCDF = std(fakeDiffCDF_Y,0,1);
% errorbar(fakeDiffCDF_X,avgFakeDiffCDF,errFakeDiffCDF)
% xlabel('Difference in Suppression')
% saveas(gcf,[image,'_DifferenceCDF.png'])
% %% test for significance in average difference
% 
% avgRealDiff = mean(realDiff);
% avgFakeDiff = sort(mean(fakeDiff,2));
% 
% 
% figure(5)
% clf
% histogram(avgFakeDiff)
% hold on
% plot([avgRealDiff,avgRealDiff],ylim)
% p = 1 - mean(avgFakeDiff < avgRealDiff);
% text(avgRealDiff,(max(ylim)/2), ['   p = ', num2str(p)])
% xlabel('Average difference in Suppression')
% saveas(gcf,[image,'_AverageDifference.png'])
% 
% 
% %% test for significance in median difference
% 
% medRealDiff = median(realDiff);
% medFakeDiff = sort(median(fakeDiff,2));
% 
% 
% figure(6)
% clf
% histogram(medFakeDiff)
% hold on
% plot([medRealDiff,medRealDiff],ylim)
% p = 1 - mean(medFakeDiff < medRealDiff);
% text(medRealDiff,(max(ylim)/2), ['   p = ', num2str(p)])
% xlabel('Median difference in suppression')
% saveas(gcf,[image,'_MedianDifference.png'])
% 
% 
% %% Plot difference accross space
% figure(7)
% clf
% hold on
% edges = linspace(0,max(dist), 100);
% 
% xmid = 0.5*(edges(1:end-1)+edges(2:end));
% [~,~,loc]=histcounts(dist,edges);
% counts = accumarray(loc,1);
% 
% %For real data
% %scatter(dist,realDiff)
% 
% meany = accumarray(loc,realDiff)./counts;
% plot(xmid,meany)
% 
% %For fake data
% %scatter(dist,fakeDiff(1,:))
% fakeYs = nan(numShuffles, length(xmid));
% 
% for i = 1:numShuffles
%     fakeYs(i,:) = accumarray(loc,fakeDiff(i,:))./counts;
% end
% 
% avgY = mean(fakeYs,1);
% errY = std(fakeYs,0,1);
% 
% errorbar(xmid,avgY,errY)
% 
% xlabel('Distance (pixels)')
% ylabel('')


%% Examine how si changes based on responsiveness
figure(8)
[~, smallSNR] = thresholdBySNR(smallGreen,smallTime,.5);
[~, bigSNR] = thresholdBySNR(bigGreen,smallTime,.5);

allSI = median(si,3);

figure(9)
clf

subplot(2,1,1)
binscatter(bigSNR(:),allSI(:))
xlabel('Big Spot SNR')
ylabel('median si')

subplot(2,1,2)
binscatter(smallSNR(:),allSI(:))
xlabel('Small Spot SNR')
ylabel('median si')



% x = smallSNR(~isnan(smallSNR));
% y = allSI(~isnan(smallSNR));
% edges = linspace(-.5,5, 100);
% [~,~,loc]=histcounts(x,edges);
% loc = loc+1;
% meany = accumarray(loc,y,[],@(x) mean(x,'omitnan'));
% meanx = accumarray(loc,x,[],@(x) mean(x,'omitnan'));
% plot(meanx,meany)

end
image