function si = calcRatioEachEpoch(smallG, bigG, smallT, bigT)
%% average accross time for each epoch
[smallGPre, smallGStim] = pullPreStim(smallG,smallT);
smallGPre = squeeze(mean(smallGPre, 3, 'omitnan'));
smallGStim = squeeze(mean(smallGStim, 3, 'omitnan'));

[bigGPre, bigGStim] = pullPreStim(bigG,bigT);
bigGPre = squeeze(mean(bigGPre, 3, 'omitnan'));
bigGStim = squeeze(mean(bigGStim, 3, 'omitnan'));

%% calculate responses (dF and dFoF) for each epoch
dF_small = smallGStim - smallGPre;
dF_big = bigGStim - bigGPre;

%% calculate SI
si = 1 - dF_big ./ dF_small;
end