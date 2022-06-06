gChannel = 2;
image = 'NoCube_00004.tif';


%% get raw image (and time image)


path = 'C:\Users\david\Desktop\UV_test';
stimChannel = 4;
start = -1;
stop = 1;
[raw,metaData] = getRaw([path filesep image]);
eImg = splitRawImageIntoEpochs2(raw(:,:,:,stimChannel), raw(:,:,:,gChannel), metaData, [start,stop]);
eResp = epochPixelResponse(eImg.green, eImg.time);

mean(eResp.stim(:)) - mean(eResp.pre(:))