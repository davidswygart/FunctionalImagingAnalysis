function [dF_SI, dFoF_SI] = calcSi_dFEpochs(smallG, bigG, smallT, bigT)
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


%% average responses accross epochs
dF_small = median(dF_small,3);
dF_big = median(dF_big,3);

dFoF_small = median(dFoF_small,3);
dFoF_big = median(dFoF_big,3);

%% calculate SI
dF_SI = (dF_small - dF_big) ./ (dF_small + dF_big);
dFoF_SI = (dFoF_small - dFoF_big) ./ (dFoF_small + dFoF_big);



end