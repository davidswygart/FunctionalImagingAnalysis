function stimInfo = StimInfo(stimTrace, numPreFrames, makePlots)
if nargin < 3
    makePlots = 0;
end

stimTrace = squeeze(stimTrace);
Deriv = diff(stimTrace);
threshold = std(Deriv);
Onset = find(Deriv > threshold) + 1;
Offset = find(Deriv < threshold*-1);  % I do not add 1 so that the offset is actually the last epoch of the stimulus ON.

Onset(diff(Onset) < 2) = []; %Remove repeats
Offset(diff(Offset) < 2) = []; %Remove repeats
numEpochs = length(Offset);
Onset = Onset(1:numEpochs); %Get rid of extra Onset if it doesn't have a corresponding offset.

if makePlots
    figure(98)
    plot(Deriv)
    hold on
    scatter(Onset, Deriv(Onset))
    scatter(Offset, Deriv(Offset))
    plot([0,length(Deriv)], [threshold, threshold])
    plot([0,length(Deriv)], [-1*threshold, -1*threshold])
    xlabel('zframes')
    ylabel('Derivative of pixel intensities')
    
    
    figure(99)
    plot(stimTrace)
    hold on
    scatter(Onset, stimTrace(Onset))
    scatter(Offset, stimTrace(Offset))
    xlabel('zframes')
    ylabel('pixelIntensity')
end

FramePerStim = Offset - Onset;
if range(FramePerStim) > 5
    error('Greater than 5 epoch difference between stim lengths')
end

FramePerStim = round(mean(FramePerStim));
FPE = round(mean(diff(Onset))); %Frames per epoch
postFrames = FPE - numPreFrames - FramePerStim;

if Offset(end) + postFrames > length(stimTrace)
   disp('remvoed last epoch because the postTime was too short')
   Onset = Onset(1:end-1);
   Offset = Offset(1:end-1);
   numEpochs = length(Offset);
end



stimInfo.numEpochs = numEpochs;
stimInfo.Onset = Onset;
stimInfo.Offset = Offset;
stimInfo.FramePerStim = FramePerStim;
stimInfo.postFrames = postFrames;
stimInfo.preFrames = numPreFrames;
stimInfo.TraceLength = FramePerStim + numPreFrames + postFrames;
end