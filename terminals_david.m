
%%
rawImage_path = 'C:/Users/david/Desktop/testImages_nonsquare/test_nonSquareImages/nonSquareImg_bi.tif';
cellName = '102721B_testImages.mat';
dataSetName = 'SpotsMultiSize_squareImg_bi';
channelOrder = ["red","green","DIC","Stim"];
%% Load Raw Image and meta-data
rawImage = getRawImage(rawImage_path, channelOrder);
%% Get threhold
zProj = mean(rawImage.img(:,:,:,1),3);
thresh = zProj > mean(zProj(:));


%% Split the image into subImages by Epoch.  Creates class instance with useful methods.

epochImg = EpochImages(rawImage, cellName, dataSetName);

%% Create image of responses accross all epochs.
epochResps = EpochResponses( )

%% plot responses for spot size 
% times = epochImg(:,:,:,5,:);
% GoR = epochImg(:,:,:,6,:);
% scatter(times(:),GoR(:))

times = epochImg(:,:,:,5,:);
G = epochImg(:,:,:,2,:);

uVals = unique(epochVal);
figure
hold on
for i = 1:length(uVals)
    vals = G(:,:,:,:,(uVals(i) == epochVal));
    t = times(:,:,:,:,(uVals(i) == epochVal));



    edges = linspace(min(t(:)),max(t(:)), 100);
    [~,~,loc]=histcounts(t(:),edges);
    meany = accumarray(loc,vals(:))./accumarray(loc,1);
    xmid = 0.5*(edges(1:end-1)+edges(2:end));

    
    plot(xmid,meany)
end

%% Motion Correction
ops.useGPU = true;
ops.kriging = true;
ops.NiterPrealign = 20;
ops.Ly = ny-1; %ignore last row, due to projector artifact

xi = ceil(.175*nx) : floor(.825*nx);
ops.Lx = numel(xi);%nx; %to cut out artifact from pockel's cell


ops = alignIterative(double(img(1:ny-1,xi,:,3)),ops); %from Suite2P ~ align using transmitted IR?
reg = rigidRegFrames(double(img(1:ny-1,:,:,1)), ops, ops.dsprealign); %from Suite2P
regIR = rigidRegFrames(double(img(1:ny-1,:,:,3)), ops, ops.dsprealign); %from Suite2P

% ops = alignIterative(reshape(epoched(1:ny-1,xi,:,:,1),ny-1, ops.Lx, []),ops); %from Suite2P ~ align using transmitted IR?
% reg = rigidRegFrames(reshape(epoched(1:ny-1,:,:,:,1),ny-1, nx, []), ops, ops.dsprealign); %from Suite2P
% regIR = rigidRegFrames(reshape(epoched(1:ny-1,xi,:,:,2),ny-1, ops.Lx, []), ops, ops.dsprealign); %from Suite2P


%% dFoF model?
ops = [];
ops.useGPU = true;
ops.fig = true;
ops.nSVDforROI = inf;
ops.diameter = 10; %!!!!
ops.yrange = 1:ny-1; ops.xrange = 1:nx;
ops.fs = frameRate;
% ops.mimg1 = mean(reshape(epoched(1:ny-1,:,mp+1:min(onTime),:),ny-1, nx, []),3);

ops.mimg1 = mean(reg,3);
% [stat,F,Fneu] = cellDetectionStandalone(reshape(epoched(1:ny-1,:,mp+1:min(onTime),:),ny-1, nx, []), ops);

[stat,F,Fneu] = cellDetectionStandalone(reg, ops);
nROIs = size(F,1);
epoched = arrayfun(@(x,y,z) cat(2,nan(nROIs,mp-x),double(F(:,z{1})), nan(nROIs,mo-y)),...
    preTime,...
    onTime,...
    ti,...
    'uniformoutput',false);
epoched = cat(3,epoched{:}) - min(F,[],'all');

bl = nanmean(epoched(:,1:mp,:),2);
dFoF = (epoched - bl) ./ (abs(bl) + eps); %this is because the noise is 0 mean...

a = splitapply(@(x) mean(x,2), shiftdim(dFoF,1), ai');

stimOn = sum(dFoF(:,mp+(1:floor(frameRate)),:).^2,2);
stimOnAvg = splitapply(@(x) mean(x,2), reshape(stimOn,nROIs,[]), ai')'; %validate this...

stimOff = sum(dFoF(:,(mp+ceil(frameRate)):end,:).^2,2);
stimOffAvg = splitapply(@(x) mean(x,2), reshape(stimOff,nROIs,[]), ai')'; %validate this...


%% Response mapping
