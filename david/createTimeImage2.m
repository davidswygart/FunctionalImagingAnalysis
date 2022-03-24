function time = createTimeImage2(sz,metadata)
nRow = sz(1);
nCol = sz(2);
nFram = sz(3);

% Does not take into account y-flyback time or reseting rows time 
inds = 1:nRow*nCol*nFram;
inds = reshape(inds,nCol, nRow, nFram);
inds = permute(inds, [2,1,3]);

% only flip rows if bidirectional scanning
if metadata.SI.hScan2D.bidirectional
    inds(2:2:end, :, :) = flip( inds(2:2:end, :, :), 2);
end

% convert indices to time values
time = inds / nRow / nCol / metadata.SI.hRoiManager.scanFrameRate;
end