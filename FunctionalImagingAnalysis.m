%% Settings
PathName = 'C:\Users\david\Desktop\GCamp6_test_images\031721B\8\'; %Path to the image
ImageName =  '031721_00008.tif'; %Image name
settings.numChannels = 3; %Number of channels in the raw image
settings.StimChannel = 3; %Which channel contains the stimulus info?
settings.AnalyzeChannel = 1; %Which channel is the flourophore (Gcamp/iGluSnFr)?
settings.preTime = 1; %How much pretime to plot (s)
settings.postTime = 1; %How much postTime to plot (s)
settings.upscaleFreq = 500; %What frequency to resample for improved stimulus alignment
settings.downscaleFreq = 50; %Downscale for easier computation
settings.respThreshold = 3; %What level of SNR must the pixel have?

%% Get Image Info
imgInfo = imfinfo([PathName ImageName]);
imgInfo = imgInfo(1);
imgInfo.SoftwareDict = char2dict(imgInfo.Software);

settings.flyback = imgInfo.SoftwareDict('SI.hScan2D.flybackTimePerFrame');
settings.pixDwell = imgInfo.SoftwareDict('SI.hScan2D.scanPixelTimeMean');
settings.frameRate = imgInfo.SoftwareDict('SI.hRoiManager.scanFrameRate');
settings.numRows = imgInfo.Height;
settings.numColumns = imgInfo.Width;
settings.XResolution = imgInfo.XResolution;
settings.YResolution = imgInfo.YResolution;


%% Load Image (Deinterleaving into stimulus and analysis images)
[image,res] = fastLoadTiff([PathName ImageName]); %Load Image
image = double(image);
image = permute(image, [2,1,3]); %Change so Image isn't rotated.  Might want to do this in fastLoadTiff()
stimImage = image(:,:,settings.StimChannel:settings.numChannels:end);
analyzeImage = image(:,:,settings.AnalyzeChannel:settings.numChannels:end);
%analyzeImage = analyzeImage - min(analyzeImage, [], 'all'); %Remove offset from anylyzeImage

%% Display The average Analysis Image (t-projection)
figure(1)
imagesc(mean(analyzeImage,3), [0,mean(analyzeImage, 'all') + .1 * std(analyzeImage, [], 'all')])
c = colorbar;
c.Label.String = 'Pixel Intensity';
title('Z projection (avgerage)')

%% Find when the stimulus turned On and Off.  Convert this to upscaled frequency.
SI = StimInfo(stimImage(end,1,:), settings); %Detect stimulus onset in orignal frames
upscaled = identifyEpochFrames(SI, settings, settings.upscaleFreq); %Find the frame numbers for PreTime, StimTime, and PostTime in upscaled frequency
downscaled = identifyEpochFrames(SI, settings, settings.downscaleFreq);

%% Calculate the delay for each pixel (in seconds)
numPixels = settings.numRows * settings.numColumns;
delay = linspace(settings.flyback, 1/settings.frameRate, numPixels); %calculate delay in seconds
%delay = round(delay*settings.upscaleFreq);%convert delay to upscaled frames
delayImage = reshape(delay, [settings.numColumns, settings.numRows]); %Reshape into the image size
delayImage = permute(delayImage, [2,1]);
delayImage(2:2:settings.numRows,:) = flip(delayImage(2:2:settings.numRows,:),2);% only needed when bidirectional scan

%% Build Image with epochs as dimension.  First use upscaled frames and then downscale them (also create SNR Image using average response during stimulus)
EpochImage = nan(settings.numRows, settings.numColumns, downscaled.numEpochFrames, SI.numEpochs);
SNRImage = nan(settings.numRows, settings.numColumns, SI.numEpochs);

for h = 1:settings.numRows
    for w = 1:settings.numColumns
        upscaledTrace = interp1(SI.frameTime, squeeze(analyzeImage(h,w,:)), upscaled.frameTime, 'previous');
        pixelDelay = delayImage(h,w);
        downscaledTrace = interp1(upscaled.frameTime+pixelDelay, upscaledTrace, downscaled.frameTime, 'linear');
        
        for e = 1:SI.numEpochs
            allFrames = downscaled.allFrames(e,:);
            preFrames = downscaled.preFrames(e,:);
            stimFrames = downscaled.stimFrames(e,:);
            postFrames = downscaled.postFrames(e,:);
            
            avgPre = mean(downscaledTrace(preFrames));
            stdPre = std(downscaledTrace(preFrames));
            avgStim = mean(downscaledTrace(stimFrames));
            SNRImage(h,w,e) = (avgStim - avgPre) / stdPre; %Calculate signal to noise ratio using preStim as noise and Stim as signal.
            
            EpochImage(h,w,:,e) = downscaledTrace(allFrames); 
        end
    end
end

%% blank out the first few rows and last row of SNR image because of bleed through.
SNRImage(1:4,:,:) = 0;
SNRImage(end,:,:) = 0;

%% Find responsive pixels (1st threshold)
avgSNRImage = mean(SNRImage, 3); %Find average SNR over all epochs for each pixel (Not splitting by spot size)

% avgSNR = mean(avgSNRImage, 'all');
% stdSNR = std(avgSNRImage, [] , 'all');
% respMask = (avgSNRImage > avgSNR+RespThreshold*stdSNR);
respMask = (avgSNRImage > settings.respThreshold); %Arbitrarily using all pixels with at least an SNR of the response Threshold

linInds = find(respMask);
[row, column] = ind2sub(size(respMask), linInds); 
matrixInds = [row, column];

%% Plot the average SNR Image and the thresholded SNR Image
figure(3)
imagesc(avgSNRImage)
c = colorbar;
c.Label.String = 'SNR';
title('SNR Image')

figure(4)
imagesc(avgSNRImage , 'AlphaData', respMask)
c = colorbar;
c.Label.String = 'SNR';
title(['First Threshold - (Threshold Mask ', num2str(settings.respThreshold), ' SNR)'])

%% Plot response trace of all responsive pixels
figure(5)
avgTrace = PlotResponseProfiles(EpochImage, matrixInds, downscaled, '1st threshold pixels');

% %% Create matched filter
% Filt = avgTrace((downscaled.numPreFrames+1):end) - mean(avgTrace(1:downscaled.numPreFrames));
% Filt = Filt / mean(Filt);
% 
% figure(6)
% plot(numPreFrames+1:length(avgTrace), Filt)
% title('matched filter')
% 
% 
% % %% Plot non threshold pixels for sanity
% % linInds = find(~respMask);
% % [nrow, ncolumn] = ind2sub(size(respMask), linInds);
% % nmatrixInds = [nrow, ncolumn];
% % avgTrace = PlotResponseProfiles(EpochImage, nmatrixInds, SI, 'NonThreshold Pixels');
% %%  Remake SNR image with Matched Filter
% SNRImage = nan(settings.numRows, settings.numChannels, SI.numEpochs);
% for h = 1:Height
%     for w = 1:Width
%         si = StimInfo(stimImage(h,w,:), numPreFrames);
%         for e = 1:SI.numEpochs
%             preFrames = si.Onset(e) - numPreFrames : si.Onset(e) - 1;
%             avgPre = mean(analyzeImage(h,w,preFrames));
%             stdPre = std(analyzeImage(h,w,preFrames));
%             
%             stimFrames = si.Onset(e):si.Onset(e)+length(Filt)-1;
%             stimData = squeeze(analyzeImage(h,w,stimFrames))';
%             filtStim = mean(Filt .* stimData);
%             
%             SNRImage(h,w,e) = (filtStim - avgPre) / stdPre;
%         end
%     end
% end
% 
% %% Cut off First few rows and last Row because of bleed through (again)
% %SNRImage = SNRImage(4:end-1,:,:); %For now just cut off first and last row of pixels
% 
% %% Find responsive pixels 
% avgSNRImage = mean(SNRImage, 3);
% 
% % avgSNR = mean(avgSNRImage, 'all');
% % stdSNR = std(avgSNRImage, [] , 'all');
% % respMask = (avgSNRImage > avgSNR+RespThreshold*stdSNR);
% respMask = (avgSNRImage > RespThreshold);
% 
% linInds = find(respMask);
% [row, column] = ind2sub(size(respMask), linInds);
% matrixInds = [row, column];
% 
% %% Plot the average SNR Image and the thresholded SNR Image
% figure(7)
% imagesc(avgSNRImage)
% c = colorbar;
% c.Label.String = 'SNR';
% title('SNR Image (from Matched Filter)')
% 
% figure(8)
% imagesc(avgSNRImage , 'AlphaData', respMask)
% c = colorbar;
% c.Label.String = 'SNR';
% title(['Second Threshold - (Threshold Mask ', num2str(settings.respThreshold), ' SNR)'])


%% determine epoch spot size (Eventually make this so that it actually reads data from Symphony.  Right now I have just been doing every other spot size)
spotSize(1:2:SI.numEpochs) = 200;
spotSize(2:2:SI.numEpochs) = 1200;
SmallInd = find(spotSize == 200);
LargeInd = find(spotSize == 1200);

%% Create Small and Large Spot SNR Images
SmallSpotSNRImage = squeeze(mean(SNRImage(:,:,SmallInd), 3));
LargeSpotSNRImage = squeeze(mean(SNRImage(:,:,LargeInd),3));


%% Find Suppression Index of all pixels above threshold
SIImage = (SmallSpotSNRImage - LargeSpotSNRImage) ./ (SmallSpotSNRImage + LargeSpotSNRImage);

figure(9)
imagesc(SIImage, 'AlphaData', respMask, [-.5,.8])
c = colorbar;
c.Label.String = 'Suppression Index (SI)';
title('SI Image (Thresholded)')

figure(10)
linInds = find(respMask);
tempVals = SIImage(linInds);
%tempVals(find(tempVals > 1)) = 1;
%tempVals(find(tempVals < -1)) = -1;
histogram(SIImage(linInds))
xlabel('Suppression Index values')
ylabel('Count')
title('SI values of thresholded pixels')

%% Plot Response over time to see beaching
EpochSNR = nan(SI.numEpochs,1);
numRespPixels = sum(respMask, 'all');
for e = 1:SI.numEpochs
    SNR_sum = sum(SNRImage(:,:,e) .* respMask, 'all');
    EpochSNR(e) = SNR_sum / numRespPixels;
end
figure(11)
plot(EpochSNR)
xlabel('epoch number')
ylabel('SNR')
title('Decrease in SNR over time')

%% Save all figures
figHandles = findall(0,'Type','figure');
savefig(figHandles, [PathName 'AllFigs'])