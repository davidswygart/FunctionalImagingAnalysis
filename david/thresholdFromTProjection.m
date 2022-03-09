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
title('time projection')

subplot(2,1,2);
imagesc(tProj, 'AlphaData', meetsThresh)
colorbar
axis image;
title('thresholded tProjection')
end