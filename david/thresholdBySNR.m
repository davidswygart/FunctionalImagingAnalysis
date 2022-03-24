function [meetsThresh, SNR] = thresholdBySNR(g,t,minSNR)
%% plot tProj of good times in the threshold channel
[pre, stim] = pullPreStim(g,t);
preAvg = squeeze(mean(pre, 3, 'omitnan'));
stimAvg = squeeze(mean(stim, 3, 'omitnan'));

difs = (stimAvg - preAvg);
difAvg = mean(difs,3);
difErr = std(difs,0,3);
SNR = difAvg ./ difErr;
%% threshold t-projection


meetsThresh = SNR > minSNR;


clf
subplot(2,1,1);
imagesc(SNR)
colorbar
axis image;
title('SNR projection')

subplot(2,1,2);
imagesc(SNR, 'AlphaData', meetsThresh)
colorbar
axis image;
title(['Meets Threshold ', num2str(minSNR)])


xlabel(['Number of pixels included = ', num2str(sum(meetsThresh(:)))])
end