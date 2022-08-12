function [dF_SI, dFoF_SI] = calcRatioPreStimEpochs(smallG, bigG, smallT, bigT)

%% average accross time for each epoch
[smallGPre, smallGStim] = pullPreStim(smallG,smallT);
smallGPre = squeeze(mean(smallGPre, 3, 'omitnan'));
smallGStim = squeeze(mean(smallGStim, 3, 'omitnan'));

[bigGPre, bigGStim] = pullPreStim(bigG,bigT);
bigGPre = squeeze(mean(bigGPre, 3, 'omitnan'));
bigGStim = squeeze(mean(bigGStim, 3, 'omitnan'));

%% average Epochs together
smallGPre = squeeze(median(smallGPre, 3, 'omitnan'));
smallGStim = squeeze(median(smallGStim, 3, 'omitnan'));

bigGPre = squeeze(median(bigGPre, 3, 'omitnan'));
bigGStim = squeeze(median(bigGStim, 3, 'omitnan'));

%% calculate SI

dF_small = smallGStim - smallGPre;
dF_big = bigGStim - bigGPre;
dF_SI = 1 - dF_big ./ dF_small;

dFoF_small = dF_small ./ smallGPre;
dFoF_big = dF_big ./ bigGPre;

dFoF_SI = 1 - dFoF_big ./ dFoF_small;



end