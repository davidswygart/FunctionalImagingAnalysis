%% Settings
PathName = 'C:\Users\david\Desktop\031721B\8\'; %Path to the image
ImageName =  '031721_00008.tif'; %Image name
preTime = 1; %How much pretime to plot (s)
FPS = 10.9; %The frame rate of the image (Hz)
numChannels = 3; %Number of channels in the raw image
StimChannel = 3; %Which channel contains the stimulus info?
AnalyzeChannel = 1; %Which channel is the flourophore (Gcamp/iGluSnFr)?
RespThreshold = 1; %What level of SNR must the pixel have?

numPreFrames = round(preTime * FPS); %How many frames to show for the preTime


%% Get image Information
info = imfinfo([PathName ImageName]);
[Height, Width] = size(imread([PathName ImageName], 1));
Zdepth = numel(info);
% Zdepth = 3000; %Temporary change to shorten processing

%% Load Image (Deinterleaving into stimulus and analysis images)
z = 1;
stimImage = zeros(Height, Width, Zdepth/numChannels);
for i=StimChannel:numChannels:Zdepth
    stimImage(:,:,z) = imread([PathName ImageName], i);
    z = z + 1;
end

z = 1;
analyzeImage = zeros(Height, Width, Zdepth/numChannels);
for i=AnalyzeChannel:numChannels:Zdepth
    analyzeImage(:,:,z) = imread([PathName ImageName], i);
    z = z + 1;
end

%analyzeImage = analyzeImage - min(analyzeImage, [], 'all'); %Remove offset from anylyzeImage

%% Display The average Analysis Image (t-projection)
figure(1)
imagesc(mean(analyzeImage,3), [0,mean(analyzeImage, 'all') + .1 * std(analyzeImage, [], 'all')])
c = colorbar;
c.Label.String = 'Pixel Intensity';
title('Z projection (avgerage)')

%% Build Image with epochs as dimension (also create SNR Image using average response during stimulus)
SI = StimInfo(stimImage(1,1,:), numPreFrames, 1); %Detect onset and offset of stimulus and plot
%SI = StimInfo(mean(AnalyzeImage(end,1:7,:), 2), numPreFrames, 1); %Temporary fix for forgotten CH4

EpochImage = nan(Height, Width, SI.TraceLength, SI.numEpochs);
EpochStim = EpochImage;
SNRImage = nan(Height, Width, SI.numEpochs);
for h = 1:Height
    for w = 1:Width
        si = StimInfo(stimImage(h,w,:), numPreFrames); % For every pixel, redetect onset and offset of stimulus.  Will be very similar to SI, but might be a frame or so shifted.
        %si = SI; %Temporary fix for forgotten CH4
        for e = 1:SI.numEpochs
            epochframes = si.Onset(e)-numPreFrames:si.Onset(e)+SI.FramePerStim+SI.postFrames-1; %Epoch Frames include preFrames, StimFrames, and Post Frames
            EpochImage(h,w,:,e) = analyzeImage(h,w,epochframes); %This is the raw analysis values for each pixel split into epochs.
            EpochStim(h,w,:,e) = stimImage(h,w,epochframes); %This image is just for sanity checks, to make sure our stimulus is lining up.
            
            preFrames = si.Onset(e) - numPreFrames : si.Onset(e) - 1; %These are the frames before the stimulus
            avgPre = mean(analyzeImage(h,w,preFrames));
            stdPre = std(analyzeImage(h,w,preFrames));
            
            stimFrames = si.Onset(e):si.Offset(e); %These are the frames during the stimulus
            avgStim = mean(analyzeImage(h,w,stimFrames));
            
            SNRImage(h,w,e) = (avgStim - avgPre) / stdPre; %Calculate signal to noise ratio using preStim as noise and Stim as signal.
        end
    end
end

%% Cut off First few rows and last Row because of bleed through
EpochImage = EpochImage(4:end-1,:,:,:);
EpochStim = EpochStim(4:end-1,:,:,:);
SNRImage = SNRImage(4:end-1,:,:); %For now just cut off first and last row of pixels



%% Find responsive pixels (1st threshold)
avgSNRImage = mean(SNRImage, 3); %Find average SNR over all epochs for each pixel (Not splitting by spot size)

% avgSNR = mean(avgSNRImage, 'all');
% stdSNR = std(avgSNRImage, [] , 'all');
% respMask = (avgSNRImage > avgSNR+RespThreshold*stdSNR);
respMask = (avgSNRImage > 5); %Arbitrarily using all pixels with at least an SNR of 5

linInds = find(respMask);
[row, column] = ind2sub(size(respMask), linInds); 
matrixInds = [row, column];

%% Plot the average SNR Image and the thresholded SNR Image
figure(2)
imagesc(avgSNRImage)
c = colorbar;
c.Label.String = 'SNR';
title('SNR Image')

figure(3)
imagesc(avgSNRImage , 'AlphaData', respMask)
c = colorbar;
c.Label.String = 'SNR';
title(['First Threshold - (Threshold Mask ', num2str(RespThreshold), ' SNR)'])

%% Create Matched Filter
figure(4)
avgTrace = PlotResponseProfiles(EpochImage, matrixInds, SI, '1st threshold pixels');
Filt = avgTrace(numPreFrames+1:end) / mean(avgTrace(numPreFrames+1:end));



% %% Plot non threshold pixels for sanity
% linInds = find(~respMask);
% [nrow, ncolumn] = ind2sub(size(respMask), linInds);
% nmatrixInds = [nrow, ncolumn];
% avgTrace = PlotResponseProfiles(EpochImage, nmatrixInds, SI, 'NonThreshold Pixels');
%%  Remake SNR image with Matched Filter
SNRImage = nan(Height, Width, SI.numEpochs);
for h = 1:Height
    for w = 1:Width
        si = StimInfo(stimImage(h,w,:), numPreFrames);
        for e = 1:SI.numEpochs
            preFrames = si.Onset(e) - numPreFrames : si.Onset(e) - 1;
            avgPre = mean(analyzeImage(h,w,preFrames));
            stdPre = std(analyzeImage(h,w,preFrames));
            
            stimFrames = si.Onset(e):si.Onset(e)+length(Filt)-1;
            stimData = squeeze(analyzeImage(h,w,stimFrames))';
            filtStim = mean(Filt .* stimData);
            
            SNRImage(h,w,e) = (filtStim - avgPre) / stdPre;
        end
    end
end

%% Cut off First few rows and last Row because of bleed through (again)
SNRImage = SNRImage(4:end-1,:,:); %For now just cut off first and last row of pixels

%% Find responsive pixels 
avgSNRImage = mean(SNRImage, 3);

% avgSNR = mean(avgSNRImage, 'all');
% stdSNR = std(avgSNRImage, [] , 'all');
% respMask = (avgSNRImage > avgSNR+RespThreshold*stdSNR);
respMask = (avgSNRImage > 5);

linInds = find(respMask);
[row, column] = ind2sub(size(respMask), linInds);
matrixInds = [row, column];

%% Plot the average SNR Image and the thresholded SNR Image
figure(5)
imagesc(avgSNRImage)
c = colorbar;
c.Label.String = 'SNR';
title('SNR Image (from Matched Filter)')

figure(6)
imagesc(avgSNRImage , 'AlphaData', respMask)
c = colorbar;
c.Label.String = 'SNR';
title(['Second Threshold - (Threshold Mask ', num2str(RespThreshold), ' stds)'])


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

figure(7)
imagesc(SIImage, 'AlphaData', respMask, [-.5,.8])
c = colorbar;
c.Label.String = 'Suppression Index (SI)';
title('SI Image (Thresholded)')

figure(8)
linInds = find(respMask);
tempVals = SIImage(linInds);
%tempVals(find(tempVals > 1)) = 1;
%tempVals(find(tempVals < -1)) = -1;
histogram(SIImage(linInds))
xlabel('Suppression Index values')
ylabel('Count')
title('SI values of thresholded pixels')

%% Plot Response over time to see beaching
EpochSNR = nan(SI.numEpochs);
for e = 1:SI.numEpochs
    EpochSNR(e) = mean(SNRImage(row, column, e), 'all');
end
figure(9)
plot(EpochSNR)
xlabel('epoch number')
ylabel('SNR')
title('Decrease in SNR over time')

%% Save all figures
figHandles = findall(0,'Type','figure');
savefig(figHandles, [PathName 'AllFigs'])