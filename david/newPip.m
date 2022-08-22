
%% get raw image (and time image)

[raw,metaData] = getRaw([path filesep image]);
%rawTime = createTimeImage2(size(raw), metaData);
%raw = raw(:,:,700:end,:);

%% Split into epochs
start = -1;
stop = 1;
eImg = splitRawImageIntoEpochs2(raw(:,:,:,stimChannel), raw(:,:,:,gChannel), metaData, [start,stop]);

%% Indicate which epochs are good
epochStart = squeeze(eImg.start);
isGoodEpoch = epochStart > minMaxTime(1) & epochStart < minMaxTime(2);

%% Show T-projection of STD
figure(1)
clf
stdG = std(eImg.green, 0, [3,4]) / abs(mean(eImg.green(:)));
imagesc(stdG)
c = colorbar;
c.Label.String = 'std';
caxis([median(min(stdG)), median(max(stdG))])
axis image
title('STD projection (only epochs) (normalized by mean)')
% saveas(gcf,[image,'_1_STD_Tprojection.png'])

% %% create rois
% rois2 = nan(100,256,12);
% plist2 = cell(12,10);
% for i = 1:12
%     p = drawpolygon()
%     rois2(:,:,i) = p.createMask();
%     plist2{i} = p;
% end


%% make threshold
threshImg = stdG > mean(stdG(:)*1.2);
figure(11)
clf
imagesc(stdG,'AlphaData',threshImg)
c = colorbar;
c.Label.String = 'std';
caxis([median(min(stdG)), median(max(stdG))])
axis image
title('threshold of STD')
% saveas(gcf,[image,'_1_STD_Tprojection.png'])

%% Get responses for each epoch (cut off bad pixels)

framResp = epochFrameResponse(eImg.green, eImg.time, threshImg);

%% Get responses for each epoch (cut off bad pixels)
% if strcmp(linescan,'x')
%     eResp = epochPixelResponse(eImg.green(:,2:end-1,:,:), eImg.time(:,2:end-1,:,:));
% else
%     %eResp = epochPixelResponse(eImg.green(1:end-1,:,:,:), eImg.time(1:end-1,:,:,:));
%     eResp = epochPixelResponse(eImg.green(16:end-1,:,:,:), eImg.time(16:end-1,:,:,:));
% end
eResp = epochPixelResponse(eImg.green, eImg.time);
% roi = mean(eImg.green(:,2:end-1,:,:),[1,2], 'omitnan');
% roiT = mean(eImg.time(:,2:end-1,:,:),[1,2], 'omitnan');

%eResp = epochPixelResponse(roi,roiT);
% imagesc(squeeze(mean(eResp.dF,3)))
% colorbar

%% Show T-projection of raw preTimes
figure(111)
clf
rawG = mean(eResp.pre(:,:,isGoodEpoch),3);
imagesc(rawG)
c = colorbar;
c.Label.String = 'image units';
caxis([median(min(rawG)), median(max(rawG))])
axis image
title('raw time projection (preTimes)')
% saveas(gcf,[image,'_1_Raw_Tprojection.png'])

%% get spot sizes for each epoch 
ne = size(eImg.green,4);
if strcmp(dataSetName, 'none')
    warning('no cell data found, assuming 2 spot sizes (200/1200)')
    spotSizes = nan(ne,1);
    spotSizes(1:2:end) = 200;
    spotSizes(2:2:end) = 1200;
else
    spotSizes = getSplitParam2(cellName,dataSetName,'curSpotSize');
    spotSizes = spotSizes(1:ne);
end

%% Plot adaption / bleaching
%framResp = squeeze(mean(eResp.dPrime, [1,2], 'omitnan'));
epochStart = squeeze(eImg.start);
uSpots = unique(spotSizes);

figure(2)
clf
title('dPrime accross epochs')
hold on
for i=1:length(uSpots)
    inds = spotSizes == uSpots(i);
    plot(epochStart(inds),framResp.dPrime(inds))
end
xlabel('epoch start time (s)')
ylabel('dPrime')

plot(xlim,[0,0], '-k')
plot(xlim,[.3,.3], '--k')
legend(num2str(uSpots),'Location','northwest')
% saveas(gcf,[image,'_2_AdaptionBleaching.png'])


figure(22)
clf
title('dFoF accross epochs')
hold on
for i=1:length(uSpots)
    inds = spotSizes == uSpots(i);
    errorbar(epochStart(inds),framResp.dFoF(inds),framResp.dFoF_sem(inds))
end
xlabel('epoch start time (s)')
ylabel('dFoF')

plot(xlim,[0,0], '--k')
%ylim([-.1,.3])
legend(num2str(uSpots),'Location','northwest')
% saveas(gcf,[image,'_2_AdaptionBleaching.png'])


%% plot sms

figure(222)
clf
spotResp = nan(length(uSpots),1);
spotErr = nan(length(uSpots),1);
for i=1:length(uSpots)
    inds = spotSizes == uSpots(i);
    spotResp(i) = mean(framResp.dFoF(inds & isGoodEpoch));
    spotErr(i) = std(framResp.dFoF(inds & isGoodEpoch));
end
errorbar(uSpots,spotResp,spotErr)
xlim([0,1250])
ylabel('average Dprime (error accross epochs)')
xlabel('spot size (um)')
hold on
plot(xlim,[0,0], '--k')
%ylim([-.1,1])
% saveas(gcf,[image,'_3_SMS.png'])
%% Choose which size for small spot

% [~,bestSizeInd] = max(spotResp(1:end-1));
% bestSize = uSpots(bestSizeInd);
[~,bestSizeInd] = min(abs(uSpots - 200));
bestSize = uSpots(bestSizeInd);
sizeStr = [num2str(round(bestSize)) ' um '];

[~, i] = max(spotResp(1:end-1));
smallInds = find((spotSizes == bestSize) & isGoodEpoch);
bigInds = find((spotSizes == uSpots(end)) & isGoodEpoch);

%drop any dangling epochs
ne = min([length(smallInds), length(bigInds)]);
smallInds = smallInds(1:ne);
bigInds = bigInds(1:ne);

%make null inds
nullInds1 = [smallInds(1:2:end); smallInds(2:2:end)];
nullInds2 = [smallInds(2:2:end); smallInds(1:2:end)];

%% show df/f for small spot
figure(3)
clf
smallDfof = eResp.dFoF(:,:,smallInds);

smallDfof = median(smallDfof,3);
imagesc(smallDfof)
c = colorbar;
c.Label.String = 'df / f';
% cmin = median(min(smallDfof));
% cmax = median(max(smallDfof));
% caxis([cmin,cmax])
caxis([-.4,4])
axis image
title([sizeStr 'dF/F'])
% saveas(gcf,[image,'_4_dFoF_small.png'])


% %% show df/f for big spot
% figure(4)
% clf
% bigDfof = eResp.dFoF(:,:,bigInds);
% bigDfof = median(bigDfof,3);
% imagesc(bigDfof)
% c = colorbar;
% c.Label.String = 'df / f';
% % cmin = min(min(smallDfof(3:end-3,3:end-3)));
% % cmax = max(max(smallDfof(3:end-3,3:end-3)));
% caxis([cmin,cmax])
% axis image
% title('Big Spot dF/F')
% 
% % saveas(gcf,[image,'_4_dFoF_big.png'])
%% show SNR (median dF / MAD) for small spot
figure(5)
clf
smallDf = eResp.dF(:,:,smallInds);
smallSNR = median(smallDf,3) ./ mad(smallDf,1,3);
imagesc(smallSNR)
c = colorbar;
c.Label.String = 'SNR (median dF / MAD)';
% cmin = median(min(smallSNR));
% cmax = median(max(smallSNR));
% caxis([cmin,cmax])
caxis([-.4,4])
axis image
title([sizeStr 'SNR'])
% saveas(gcf,[image,'_5_SNR_small.png'])

% %% show SNR (median dF / MAD) for big spot
% figure(6)
% clf
% bigDf = eResp.dF(:,:,bigInds);
% bigSNR = median(bigDf,3) ./ mad(bigDf,1,3);
% imagesc(bigSNR)
% c = colorbar;
% c.Label.String = 'SNR (median dF / MAD)';
% 
% % cmin = mean(min(bigSNR(4:end-4,4:end-4)));
% % cmax = mean(max(bigSNR(4:end-4,4:end-4)));
% caxis([cmin,cmax])
% axis image
% title('Big Spot SNR')
% 
% saveas(gcf,[image,'_6_SNR_big.png'])

%% show average dPrime for small spot
figure(6)
clf
smallDprime = eResp.dPrime(:,:,smallInds);
smallDprime = median(smallDprime,3,'omitnan');
imagesc(smallDprime)
c = colorbar;
c.Label.String = 'median dPrime';
% cmin = median(min(smallSNR));
% cmax = median(max(smallSNR));
% caxis([cmin,cmax])
caxis([-.2,2])
axis image
title([sizeStr 'Dprime'])
% saveas(gcf,[image,'_6_Dprime_small.png'])

%% get response and null values
smallDf = eResp.dF(:,:,smallInds);
bigDf = eResp.dF(:,:,bigInds);

nullDf1 = eResp.dF(:,:,nullInds1);
nullDf2 = eResp.dF(:,:,nullInds2);

%% plot si
[medSI,allSI] = calcAndPlotSi(smallDf,bigDf, 1);
%medSI = (1 - medSI)*100;
figure(7)
clf
imagesc(medSI)
c = colorbar;
c.Label.String = '(1 - Big/Small)*100';
caxis([-0,150])
axis image
title('median suppression')
%saveas(gcf,[image,'_7_si_real.png'])

figure(77)
clf
imagesc(medSI,'AlphaData',threshImg)
c = colorbar;
c.Label.String = '(1 - Big/Small)*100';
caxis([-0,150])
axis image
title('median suppression')
%saveas(gcf,[image,'_7_si_real.png'])
%%

figure(777)
clf
imagesc(imgaussfilt(medSI,1),'AlphaData',threshImg)
c = colorbar;
c.Label.String = '(1 - Big/Small)*100';
caxis([-0,150])
axis image
title('median suppression')
%saveas(gcf,[image,'_7_si_real.png'])

%% plot null si
% figure(15)
% clf
% [siN,siNall] = calcAndPlotSi(nullDf1,nullDf2, 1);
% %saveas(gcf,[image,'_8_si_null.png'])
% % % plot null si error
% % figure(16)
% % clf
% siN_mad = mad(nullDf2 ./ nullDf1, 1,3);