function avgTrace = PlotResponseProfiles(EpochImage, matrixInds, SI, Title)
[~,~,z,~] = size(EpochImage);
numPixels = length(matrixInds);
traces = nan(numPixels, z);
for p = 1:numPixels
    Pixel = EpochImage(matrixInds(p,1), matrixInds(p,2),:,:);
    trace = squeeze(mean(Pixel, 4));
    plot(trace)
    hold on
    traces(p,:) = trace;
end

avgTrace = squeeze(mean(traces,1));

hold on
plot(avgTrace, 'LineWidth', 6)
hold on
plot([SI.preFrames+1, SI.preFrames+SI.FramePerStim+1], [0,0], 'y',  'LineWidth', 6)
ylabel('Pixel Values')
xlabel('frames')
title(Title)