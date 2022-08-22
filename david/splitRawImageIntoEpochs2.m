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


frameStart = on(:,3) - nPreFram;
frameEnd = on(:,3) + nPostFram;

if frameStart(1) < 0
    warning('The first image epoch occurs earlier than the pretime length.  Removing this epoch.')
    frameEnd(frameStart < 0) = [];
    frameStart(frameStart < 0) = [];
end

if frameEnd(end) > size(green,3)
    warning('The last image epoch postime extends beyond the imaging time.  Removing this epoch.')
    frameStart(frameEnd > size(green,3)) = [];
    frameEnd(frameEnd > size(green,3)) = [];
end

ne = length(frameStart);

eImg.green = nan(nRow, nCol, nTotalFram, ne);
eImg.time = eImg.green;
eImg.stim = eImg.green;
eImg.start = nan(1,1,1,ne);

for en = 1:ne       
    r = on(en,1);
    c = on(en,2);
    f = on(en,3); 
    framInds = frameStart(en):frameEnd(en);
    
    eImg.green(:,:,:,en) = green(:,:,framInds);
    eImg.stim(:,:,:,en) = stim(:,:,framInds);
    
    time = rawTime(:,:,framInds);
    eImg.start(en) = time(1);
    eImg.time(:,:,:,en) = time - rawTime(r,c,f);
end

%% plot avg stim to verify the epochs were split correctly
figure(1111)
clf
plot(squeeze(mean(eImg.time, [1,2,4])), squeeze(mean(eImg.stim, [1,2,4])))
%plot(squeeze(mean(epochs.time, [1,2,4])), squeeze(mean(epochs.green, [1,2,4])))
end
