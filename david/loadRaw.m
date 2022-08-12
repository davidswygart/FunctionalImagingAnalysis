path = 'C:\Users\david\Desktop\2p_copies';

%% load the raw image and metadata
[tempImg, raw.res, raw.md, ~,raw.position] = fastLoadTiff([path filesep image]); %not using timestamps
tempImg = permute(tempImg, [2, 1, 4, 3]);

%% organize the raw data and metadata into a struct (raw)
raw.green = tempImg(:,:,:,gChannel);
raw.stim = tempImg(:,:,:,stimChannel);
raw.thresh = tempImg(:,:,:,threshChannel);

% parse metatdata into struct to be more readable
s = parseScanImageMetaData(raw.md);
raw.SI = s.SI;

% set some commonly used metatdata values to top of struct for convenience
[raw.nRow, raw.nCol, raw.nFram] = size(raw.green);
raw.framRat = raw.SI.hRoiManager.scanFrameRate;
raw.isBi = raw.SI.hScan2D.bidirectional;

%% make a time image
raw.time = createTimeImage(raw.nRow, raw.nCol, raw.nFram, raw.framRat,raw.isBi);
raw.framT = squeeze(raw.time(1,1,:));

%% subtract off the mode and get rid of 0 values
histogram(raw.green(:))

s = sort(raw.green(:));
sub = s(round(length(s) * .001)); % find the dimmest 0.1% pixel to subtract off of the image
raw.green = raw.green - sub;
raw.green(raw.green < 0) = 0;

%% Set the end pixels to 0 (artifact) -> applied to both green and threshold channels
if linescan == 'y'
    raw.green([1,2,3,end],:,:) = 0;
    raw.thresh([1,2,3,end],:,:) = 0;
else
    raw.green(:,[1,end],:) = 0;
    raw.thresh(:,[1,end],:) = 0;
end

%%
%raw.green = randn(size(raw.green));

