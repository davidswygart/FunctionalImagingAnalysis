function avgTrace = PlotResponseProfiles(EpochImage, matrixInds, frames, Title)
[~,~,z,~] = size(EpochImage);
numPixels = length(matrixInds);
traces = nan(numPixels, z);
for p = 1:numPixels
    Pixel = EpochImage(matrixInds(p,1), matrixInds(p,2),:,:);
    trace = squeeze(mean(Pixel, 4));
    plot(frames.epochFrameTime, trace)
    hold on
    traces(p,:) = trace;
end

avgTrace = squeeze(mean(traces,1));

hold on
plot(frames.epochFrameTime, avgTrace, 'LineWidth', 6)
hold on

onTime = frames.epochFrameTime(frames.numPreFrames+1);
offTime = frames.epochFrameTime(frames.numPreFrames+frames.numStimFrames);
plot([onTime, offTime], [0,0], 'y',  'LineWidth', 8) %plot a bar to indicate the stimulus
ylabel('Pixel Values')
xlabel('frames')
title(Title)