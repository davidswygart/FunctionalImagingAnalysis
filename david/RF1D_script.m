path = 'D:\Images\032422B';
image = 'rf1d_00006.tif';

cellName = 'none';
dataSetName = 'none';

gChannel = 2;
stimChannel = 4;
threshChannel = 2;
minMaxTime = [0 inf];
linescan = 'x';

start = -1;
stop = 1;

numPos = 7;
barSep = 60;
repeats = 3;


%% get raw image (and time image)
[raw,metaData] = getRaw([path filesep image]);
%rawTime = createTimeImage2(size(raw), metaData);


%% Split into epochs
figure(111)
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

%% get RF1D location

spotSizes = genRF1Dparams(numPos,barSep,repeats);
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
spotErr = nan(length(uSpots),1);
for i=1:length(uSpots)
    inds = spotSizes == uSpots(i);
    spotResp(i) = mean(framResp(inds & isGoodEpoch));
    spotErr(i) = std(framResp(inds & isGoodEpoch));
end
[~,bestSizeInd] = max(spotResp(1:end-1));
bestSize = uSpots(bestSizeInd);

[~, i] = max(spotResp(1:end-1));
smallInds = find((spotSizes == bestSize) & isGoodEpoch);

%% Plot responses by RF1D position
figure(30)
clf
errorbar(uSpots,spotResp,spotErr)
xlabel('position')
ylabel('dPrime')

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
cmin = 0;
cmax = 2;

caxis([cmin,cmax])
axis image
title('Small Spot dF/F')

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
cmin = 0;
cmax = 2;
caxis([cmin,cmax])
axis image
title('Small Spot SNR')

%%
figure(30)