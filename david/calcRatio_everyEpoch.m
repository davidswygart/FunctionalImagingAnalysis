function [dF_SI, dFoF_SI] = calcRatio_everyEpoch(smallG, bigG, smallT, bigT)
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

dFoF_small = dF_small ./ smallGPre;
dFoF_big = dF_big ./ bigGPre;

%% calculate SI
dF_SI = 1 - dF_big ./ dF_small;
dFoF_SI = 1 - dFoF_big ./ dFoF_small;

%% average SI accross epochs
dF_SI = median(dF_SI,3);
dFoF_SI = median(dFoF_SI,3);
end