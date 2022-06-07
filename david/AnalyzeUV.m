function res = AnalyzeUV()
path = 'D:\Images\060622B_PMT_filterTest';
image = 'doubleGreenBandAndDavidUV_00001';

stimChannel = 4;
start = -1;
stop = 2;
%% get raw image (and time image)
[raw,metaData] = getRaw([path filesep image '.tif']);
eImgCh1 = splitRawImageIntoEpochs2(raw(:,:,:,stimChannel), raw(:,:,:,1), metaData, [start,stop]);
eImgCh2 = splitRawImageIntoEpochs2(raw(:,:,:,stimChannel), raw(:,:,:,2), metaData, [start,stop]);

eRespCh1 = epochPixelResponse(eImgCh1.green, eImgCh1.time);
eRespCh2 = epochPixelResponse(eImgCh2.green, eImgCh2.time);

avgCh1 = squeeze(mean(eRespCh1.stim, [1,2]) - mean(eRespCh1.pre, [1,2]));
avgCh2 = squeeze(mean(eRespCh2.stim, [1,2]) - mean(eRespCh2.pre, [1,2]));

res = struct;
res.avg = [mean(avgCh1) mean(avgCh2)];
res.std = [std(avgCh1) std(avgCh2)];
end