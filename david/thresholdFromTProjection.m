function [meetsThresh, tProj] = thresholdFromTProjection(img, p)
%% plot tProj of good times in the threshold channel
tProj = mean(img,3,'omitnan');



%% threshold t-projection
s = sort(tProj(:));
minV = s(round(length(s)*(1-p)));

meetsThresh = tProj > minV;


clf
subplot(2,1,1);
imagesc(tProj)
colorbar
axis image;
title('Intensity Projection through time')

subplot(2,1,2);
imagesc(tProj, 'AlphaData', meetsThresh)
colorbar
axis image;
title(['Meets Threshold ', num2str(p)])


xlabel(['Number of pixels included = ', num2str(sum(meetsThresh(:)))])
end