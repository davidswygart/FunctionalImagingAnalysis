function [dF_SI, dFoF_SI] = calcRatioRaw(smallG, bigG, smallT, bigT)

[smallGPre, smallGStim] = pullPreStim(smallG,smallT);
smallGPre = squeeze(median(smallGPre, [3,4], 'omitnan'));
smallGStim = squeeze(median(smallGStim, [3,4], 'omitnan'));

[bigGPre, bigGStim] = pullPreStim(bigG,bigT);
bigGPre = squeeze(median(bigGPre, [3,4], 'omitnan'));
bigGStim = squeeze(median(bigGStim, [3,4], 'omitnan'));

dF_small = smallGStim - smallGPre;
dF_big = bigGStim - bigGPre;
dF_SI = 1 - dF_big ./ dF_small;

dFoF_small = dF_small ./ smallGPre;
dFoF_big = dF_big ./ bigGPre;

dFoF_SI = 1 - dFoF_big ./ dFoF_small;
end