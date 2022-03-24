
%% get raw image (and time image)
[raw,metaData] = getRaw(image);
%rawTime = createTimeImage2(size(raw), metaData);


%% Split into epochs
start = -1;
stop = 2;
eImg = splitRawImageIntoEpochs2(raw(:,:,:,stimChannel), raw(:,:,:,gChannel), metaData, [start,stop]);
%% Get responses for each epoch (cut off bad pixels)
if strcmp(linescan,'x')
    eResp = epochPixelResponse(eImg.green(:,2:end-1,:,:), eImg.time(:,2:end-1,:,:));
else
    %eResp = epochPixelResponse(eImg.green(1:end-1,:,:,:), eImg.time(1:end-1,:,:,:));
    eResp = epochPixelResponse(eImg.green(16:end-1,:,:,:), eImg.time(16:end-1,:,:,:));
end


%% Indicate which epochs are good
epochStart = squeeze(eImg.start);
isGoodEpoch = epochStart > minMaxTime(1) & epochStart < minMaxTime(2);
%% Show T-projection of raw preTimes
figure(1)
clf
rawG = mean(eResp.pre(:,:,isGoodEpoch),3);
imagesc(rawG)
c = colorbar;
c.Label.String = 'image units';
caxis([median(min(rawG)), median(max(rawG))])
axis image
title('raw time projection (preTimes)')
saveas(gcf,[image,'_1_Raw_Tprojection.png'])

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
framResp = squeeze(mean(eResp.dPrime, [1,2], 'omitnan'));
epochStart = squeeze(eImg.start);
uSpots = unique(spotSizes);

figure(2)
clf
hold on
for i=1:length(uSpots)
    inds = spotSizes == uSpots(i);
    plot(epochStart(inds),framResp(inds))
end
xlabel('epoch start time (s)')
ylabel('average dPrime')
legend(num2str(uSpots),'Location','northwest')
saveas(gcf,[image,'_2_AdaptionBleaching.png'])




%% Determine best spot size (median of mean dFoF for good epochs)from adaptation plot
spotResp = nan(length(uSpots),1);
for i=1:length(uSpots)
    inds = spotSizes == uSpots(i);
    spotResp(i) = median(framResp(inds & isGoodEpoch));
end
[~,bestSizeInd] = max(spotResp(1:end-1));
bestSize = uSpots(bestSizeInd);

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
cmin = median(min(smallDfof));
cmax = median(max(smallDfof));
caxis([cmin,cmax])
axis image
title('Small Spot dF/F')
saveas(gcf,[image,'_3_dFoF_small.png'])
%% show df/f for big spot
figure(4)
clf
bigDfof = eResp.dFoF(:,:,bigInds);
bigDfof = median(bigDfof,3);
imagesc(bigDfof)
c = colorbar;
c.Label.String = 'df / f';
% cmin = min(min(smallDfof(3:end-3,3:end-3)));
% cmax = max(max(smallDfof(3:end-3,3:end-3)));
caxis([cmin,cmax])
axis image
title('Big Spot dF/F')

saveas(gcf,[image,'_4_dFoF_big.png'])
%% show SNR (median dF / MAD) for small spot
figure(5)
clf
smallDf = eResp.dF(:,:,smallInds);
smallSNR = median(smallDf,3) ./ mad(smallDf,1,3);
imagesc(smallSNR)
c = colorbar;
c.Label.String = 'SNR (median dF / MAD)';
cmin = median(min(smallSNR));
cmax = median(max(smallSNR));
caxis([cmin,cmax])
axis image
title('Small Spot SNR')
saveas(gcf,[image,'_5_SNR_small.png'])
%% show SNR (median dF / MAD) for big spot
figure(6)
clf
bigDf = eResp.dF(:,:,bigInds);
bigSNR = median(bigDf,3) ./ mad(bigDf,1,3);
imagesc(bigSNR)
c = colorbar;
c.Label.String = 'SNR (median dF / MAD)';

% cmin = mean(min(bigSNR(4:end-4,4:end-4)));
% cmax = mean(max(bigSNR(4:end-4,4:end-4)));
caxis([cmin,cmax])
axis image
title('Big Spot SNR')

saveas(gcf,[image,'_6_SNR_big.png'])
