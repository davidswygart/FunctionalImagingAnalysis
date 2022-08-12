function pos = genRF1Dparams(numPos,barSep,repeats)
    pos = zeros(numPos,1);
    steps = barSep:barSep:(barSep*(numPos-1)/2);
    pos(2:2:end) = steps*-1;
    pos(3:2:end) = steps;
    pos = repmat(pos,repeats,1);
end