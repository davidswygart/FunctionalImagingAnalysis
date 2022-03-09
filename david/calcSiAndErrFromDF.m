function [si, err] = calcSiAndErrFromDF(smallG, bigG, smallT, bigT)
%% average accross time for each epoch
[smallGPre, smallGStim] = pullPreStim(smallG,smallT);
smallGPre = squeeze(mean(smallGPre, 3, 'omitnan'));
smallGStim = squeeze(mean(smallGStim, 3, 'omitnan'));

[bigGPre, bigGStim] = pullPreStim(bigG,bigT);
bigGPre = squeeze(mean(bigGPre, 3, 'omitnan'));
bigGStim = squeeze(mean(bigGStim, 3, 'omitnan'));

%% calculate responses (dF) for each epoch
dF_small = smallGStim - smallGPre;
dF_big = bigGStim - bigGPre;

%% average responses accross epochs
dF_small_avg = mean(dF_small,3);
dF_big_avg = mean(dF_big,3);

%% calculate SI
si = (dF_small_avg - dF_big_avg) ./ (dF_small_avg + dF_big_avg);

%% error of dF
dF_small_err = std(dF_small,0,3);
dF_big_err = std(dF_big,0,3);

%% error of SI
num = dF_small_avg - dF_big_avg;
den = dF_small_avg + dF_big_avg;

num_err = sqrt(dF_small_err.^2 + dF_big_err.^2);
den_err = num_err;
err = abs(si) .* sqrt((num_err ./ num).^2 + (den_err ./ den).^2);
end