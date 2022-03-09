function eImg = splitRawImageIntoEpochs(raw)
%% find epoch start from stimulus channel
[on, off] = detectStimInds(raw.stim, raw.isBi);
numOn = size(on,1);

%%  create epoch image (nrows, nColumns, nTframes, nEpochs) EpochImage & EpochTimes
nPreFram = ceil(1 * raw.framRat);
nPostFram = ceil(2 * raw.framRat);
nTotalFram = nPreFram+nPostFram+1;

eImg = struct;
eImg.green = nan(raw.nRow, raw.nCol, nTotalFram, numOn);
eImg.time = eImg.green;
eImg.stim = eImg.green;
eImg.start = nan(1,1,1,numOn);
for en = 1:numOn       
    r = on(en,1);
    c = on(en,2);
    f = on(en,3); 
    framInds = (f-nPreFram):(f+nPostFram);
    
    eImg.green(:,:,:,en) = raw.green(:,:,framInds);
    eImg.stim(:,:,:,en) = raw.stim(:,:,framInds);
    
    time = raw.time(:,:,framInds);
    eImg.start(en) = time(1);
    eImg.time(:,:,:,en) = time - raw.time(r,c,f);
end

%% plot avg stim to verify the epochs were split correctly
clf
plot(squeeze(mean(eImg.time, [1,2,4])), squeeze(mean(eImg.stim, [1,2,4])))
%plot(squeeze(mean(epochs.time, [1,2,4])), squeeze(mean(epochs.green, [1,2,4])))
end
