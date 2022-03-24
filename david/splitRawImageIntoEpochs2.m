function eImg = splitRawImageIntoEpochs2(stim, green, metaData, range)
%% find epoch start from stimulus channel
[on, off] = detectStimInds(stim, metaData.SI.hScan2D.bidirectional);
numOn = size(on,1);

%% make a time image
rawTime = createTimeImage2(size(stim), metaData);

%%  create epoch image (nrows, nColumns, nTframes, nEpochs) EpochImage & EpochTimes
framRat = metaData.SI.hRoiManager.scanFrameRate;

nPreFram = ceil(-1*range(1) * framRat);
nPostFram = ceil(range(2) * framRat);
nTotalFram = nPreFram+nPostFram+1;

[nRow,nCol,~] = size(green);

eImg = struct;
eImg.green = nan(nRow, nCol, nTotalFram, numOn);
eImg.time = eImg.green;
eImg.stim = eImg.green;
eImg.start = nan(1,1,1,numOn);
for en = 1:numOn       
    r = on(en,1);
    c = on(en,2);
    f = on(en,3); 
    framInds = (f-nPreFram):(f+nPostFram);
    
    eImg.green(:,:,:,en) = green(:,:,framInds);
    eImg.stim(:,:,:,en) = stim(:,:,framInds);
    
    time = rawTime(:,:,framInds);
    eImg.start(en) = time(1);
    eImg.time(:,:,:,en) = time - rawTime(r,c,f);
end

%% plot avg stim to verify the epochs were split correctly
clf
plot(squeeze(mean(eImg.time, [1,2,4])), squeeze(mean(eImg.stim, [1,2,4])))
%plot(squeeze(mean(epochs.time, [1,2,4])), squeeze(mean(epochs.green, [1,2,4])))
end
