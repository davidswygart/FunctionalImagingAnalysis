function time = createTimeImage(nRow,nCol,nFram,framRat,isBi)
% Does not take into account y-flyback time or reseting rows time 
inds = 1:nRow*nCol*nFram;
inds = reshape(inds,nCol, nRow, nFram);
inds = permute(inds, [2,1,3]);

% only flip rows if bidirectional scanning
if isBi
    inds(2:2:end, :, :) = flip( inds(2:2:end, :, :), 2);
end

% convert indices to time values
time = inds / nRow / nCol / framRat;
end