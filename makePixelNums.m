function pixNums = makePixelNums(sz, bidirectional)
    lin = 1:int64(sz(1)*sz(2)*sz(3));

    pixNums = reshape(lin,[sz(2),sz(1),sz(3)]);
    pixNums = permute(pixNums,[2,1,3]);

    if(bidirectional)
        pixNums(2:2:end,:,:) = flip(pixNums(2:2:end,:,:),2);
    end
end