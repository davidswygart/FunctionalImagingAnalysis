
%
% rawImage_path = 'C:/Users/david/Desktop/102121B/102121Bc2/Raw/102121Bc2_SMS.tif';
% cellName = '102121Bc2';
% dataSetName = 'SpotsMultiSize_functionalImaging';
% channelOrder = ["green","red","DIC","Stim"];

% rawImage_path = 'C:/Users/david/Desktop/testImages_nonsquare/test_nonSquareImages/nonSquareImg_bi.tif';
% cellName = '102721B_testImages.mat';
% dataSetName = 'SpotsMultiSize_squareImg_bi';
% channelOrder = ["red","green","DIC","Stim"];

rawImage_path = 'C:/Users/david/Desktop/021022B/cell1/SMS.tif';
cellName = '021022Bc1.mat';
dataSetName = 'SpotsMultiSize';
channelOrder = ["red","green","DIC","Stim"];

%% Load raw image class
rImg = RawImage(rawImage_path);



%%
%rImg.detectStimLines(4)
rImg.createEpochImage(cellName, dataSetName, 4)

%% Load Raw Image and meta-data and reset 0
rawImage = getRawImage(rawImage_path, channelOrder);
%set new 0 for green and red based on minimum mode;
%ms = (mode(rawImage.img(:,:,:,1:2),3)); %mode for every pixel for red/green
%minMod = min(ms,[],[1 2]); %min pixel mode for red/green
% rawImage.img(:,:,:,1) = rawImage.img(:,:,:,1) - minMod(1); %subtrace mode for ch1
% rawImage.img(:,:,:,2) = rawImage.img(:,:,:,2) - minMod(2); %subtrace mode for ch2
rawImage.img(:,:,:,1:2) = rawImage.img(:,:,:,1:2) + 32768; %undo the calibration from scanimage
rawImage.img = uint64(rawImage.img);% convert to unsigned integer to cut off negative values
%% Split the image into subImages by Epoch.  Creates class instance with useful methods.
epochImg = EpochImages(rawImage, cellName, dataSetName);

%% Get threhold
figure
zProj = mean(rawImage.img(:,:,:,2),3);
thresh = zProj > mean(zProj(:));
%thresh = zProj > 300;
scnan(zProj, thresh);



%% Analysis
[SNR, dFoF, avgStim, avgPre, dF, preVals, stimVals] = epochImg.epochResponses();
%[SNR, dFoF, avgStim, avgPre] = epochImg.epochResponses('channel', "green", 'ratioChannel', "red");



%% Average all epochs together before calclulating dF/SNR
spotsizeInd = 5;

[eIdentity, spotSizes] = epochImg.epochsByParam('curSpotSize');
smallEpochs = find(eInd == spotsizeInd); %only 100 um
%plotEpochs = find(eIdentity); %all spot sizes

avgAvgPre = squeeze(mean(preVals(:,:,:,smallEpochs),[3 4],'omitnan')); %average accross time and epochs
avgAvgStim = squeeze(mean(stimVals(:,:,:,smallEpochs),[3 4],'omitnan'));


avgAvgDf = avgAvgStim - avgAvgPre;

figure('Name','avgAvgDf')
imagesc(avgAvgDf)
colorbar

avgAvgSTDPre = squeeze(std(preVals(:,:,:,smallEpochs), 0,[3 4],'omitnan')); %std accross time and epochs (weight 0)
figure('Name','avgAvgSTD')
imagesc(avgAvgSTDPre)
colorbar

avgAvgSNR = avgAvgDf ./ avgAvgSTDPre;

figure('Name','avgAvgSNR_smallSpot')
imagesc(avgAvgSNR)
colorbar


%% Use above to calculate SI
SI = (sSNR - lSNR) ./ (sSNR + lSNR);

figure('Name','avgAvgSI')
scnan(SI, SNRThresh)
caxis([-1 3])
colorbar

figure('Name','histogramSI')
edges = [-inf linspace(-1,2,100) inf];
h = histogram(SI(SNRThresh),edges)
%% threhold by SNR
SNRThresh = avgAvgSNR > .2;
scnan(avgAvgSNR, SNRThresh)
colorbar
caxis([-2 2])

figure
edges = [-inf linspace(-3,3,100) inf];
histogram(SI(SNRThresh))
%% bleaching of accross epochs
data = avgStim;
spotsizeInd = 2;

[eIdentity, spotSizes] = epochImg.epochsByParam('curSpotSize');
plotEpochs = find(eInd == spotsizeInd); %only pick certain spot size
%plotEpochs = find(eIdentity); %all spot sizes

data = mean(data,[1 2]);%average data accross xy
%data = data(28,83,:);%for specific pixel (switch XY)

x = plotEpochs';
y = squeeze(data(plotEpochs));
scatter(x,y)
f = polyfit(x,y,1); %fit line

%% mean SNR accross all epochs
mSNR = squeeze(mean(SNR,3));
%mSNR(~thresh) = 0;
imagesc(mSNR)
colorbar

%% mean dFoF accross all epochs
mdFoF = squeeze(mean(dFoF,3));
%mdFoF(~thresh) = 0;
imagesc(mdFoF)
colorbar

%% mean of above values for spot size 100
data = SNR;

[eInd, uniqueVals] = epochImg.epochsByParam('curSpotSize');
es = eInd == 2;

data = data(:,:,es);

% pTiles = prctile(data,[20 80], 3);
% pinds = (data > pTiles(:,:,1))  &   (data < pTiles(:,:,2));

stdData = std(data,0, 3);
mData = mean(data, 3);
mult = 0.00000001;
scutoff = (stdData * mult) - mData;
lcutoff = (stdData * mult) + mData;
pinds = (data > scutoff)  &   (data < lcutoff);

tempData = data;
tempData(~pinds) = nan;

mSNR = squeeze(mean(tempData,3, 'omitnan'));


%mSNR(~thresh) = 0;
imagesc(mSNR)
colorbar

figure
plotPoints = [10,99; 32,57; 61,75; 100,101];
hold on
for i = 1:4
    subplot(4,1,i)
    hold on
    xind = plotPoints(i,1);
    yind = plotPoints(i,2);
    zinds = pinds(xind,yind,:);
    plot(squeeze(data(xind,yind,zinds)));
    plot([0 sum(zinds)], [0 0])
end
xlabel('epoch')
subplot(4,1,2)
ylabel('SNR')
hold off



%% find infinite
big = (SNR == inf);
[r,c,v] = find(big); % No positive infinite

small = (SNR == -inf);
[x,y,z] = ind2sub(size(SNR), find(small));
ind = [x,y,z];

badPix = zeros(size(SNR));
badPix(small) = 1;
badProj = squeeze(sum(badPix,3));
imagesc(badProj);

badProj(~thresh) = 0;
imagesc(badProj);

%% Find nan Pix
n = isnan(SNR);
[x,y,z] = ind2sub(size(SNR), find(n));
ind = [x,y,z];

badPix = zeros(size(SNR));
badPix(n) = 1;
badProj = squeeze(sum(badPix,3));
imagesc(badProj);

badProj(~thresh) = 0;
imagesc(badProj);

%% Get ROI
%ROIS - 2D image where # ROI num
%params = {'paramNameAndValue',{'curSpotSize',110.668}};
%params = {'paramNameAndValue',{'curSpotSize',110.68}};


params = {'curSpotSize'};
r = ROIs(@epochImg.plotSplitByParam, params);
%r.addROI();

%% plot data against time
epochImg.plotRaw(); %big response to green
epochImg.plotRaw('ChannelY', 'red'); %small response to green
epochImg.plotSplitByParam('curSpotSize'); %size 110 um best
epochImg.plotSplitByParam('curSpotSize', 'xyMask',thresh); %Similar responses but just larger when adding a mask


%epochImg.plotOverTrials('paramNameAndValue',{'curSpotSize',110.668});%some bleaching
%%
% %% Motion Correction
% ops.useGPU = true;
% ops.kriging = true;
% ops.NiterPrealign = 20;
% ops.Ly = ny-1; %ignore last row, due to projector artifact
% 
% xi = ceil(.175*nx) : floor(.825*nx);
% ops.Lx = numel(xi);%nx; %to cut out artifact from pockel's cell
% 
% 
% ops = alignIterative(double(img(1:ny-1,xi,:,3)),ops); %from Suite2P ~ align using transmitted IR?
% reg = rigidRegFrames(double(img(1:ny-1,:,:,1)), ops, ops.dsprealign); %from Suite2P
% regIR = rigidRegFrames(double(img(1:ny-1,:,:,3)), ops, ops.dsprealign); %from Suite2P
% 
% % ops = alignIterative(reshape(epoched(1:ny-1,xi,:,:,1),ny-1, ops.Lx, []),ops); %from Suite2P ~ align using transmitted IR?
% % reg = rigidRegFrames(reshape(epoched(1:ny-1,:,:,:,1),ny-1, nx, []), ops, ops.dsprealign); %from Suite2P
% % regIR = rigidRegFrames(reshape(epoched(1:ny-1,xi,:,:,2),ny-1, ops.Lx, []), ops, ops.dsprealign); %from Suite2P
% 
% 
% %% dFoF model?
% ops = [];
% ops.useGPU = true;
% ops.fig = true;
% ops.nSVDforROI = inf;
% ops.diameter = 10; %!!!!
% ops.yrange = 1:ny-1; ops.xrange = 1:nx;
% ops.fs = frameRate;
% % ops.mimg1 = mean(reshape(epoched(1:ny-1,:,mp+1:min(onTime),:),ny-1, nx, []),3);
% 
% ops.mimg1 = mean(reg,3);
% % [stat,F,Fneu] = cellDetectionStandalone(reshape(epoched(1:ny-1,:,mp+1:min(onTime),:),ny-1, nx, []), ops);
% 
% [stat,F,Fneu] = cellDetectionStandalone(reg, ops);
% nROIs = size(F,1);
% epoched = arrayfun(@(x,y,z) cat(2,nan(nROIs,mp-x),double(F(:,z{1})), nan(nROIs,mo-y)),...
%     preTime,...
%     onTime,...
%     ti,...
%     'uniformoutput',false);
% epoched = cat(3,epoched{:}) - min(F,[],'all');
% 
% bl = nanmean(epoched(:,1:mp,:),2);
% dFoF = (epoched - bl) ./ (abs(bl) + eps); %this is because the noise is 0 mean...
% 
% a = splitapply(@(x) mean(x,2), shiftdim(dFoF,1), ai');
% 
% stimOn = sum(dFoF(:,mp+(1:floor(frameRate)),:).^2,2);
% stimOnAvg = splitapply(@(x) mean(x,2), reshape(stimOn,nROIs,[]), ai')'; %validate this...
% 
% stimOff = sum(dFoF(:,(mp+ceil(frameRate)):end,:).^2,2);
% stimOffAvg = splitapply(@(x) mean(x,2), reshape(stimOff,nROIs,[]), ai')'; %validate this...
