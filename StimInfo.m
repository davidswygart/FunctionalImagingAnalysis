function outputStruct = StimInfo(stimTrace, settings)
%% Find Onset and offset of stimulus (in z indices)
stimTrace = squeeze(stimTrace);
Deriv = diff(stimTrace);
threshold = std(Deriv);
Onset = find(Deriv > threshold) + 2; % I add 1 because everything is shifted 1 by diff(), I add another 1 because the frame must be on before the stimulus is actually shown during the flyback
Offset = find(Deriv < threshold*-1) + 2;

Onset(diff(Onset) < 2) = []; %Remove repeats
Offset(diff(Offset) < 2) = []; %Remove repeats
Onset = Onset(1:length(Offset)); %Get rid of extra Onset if it doesn't have a corresponding offset.

%% Find average length of the stimulus (easier to just go with average so that traces are consistent, but we are giving up a little fidelity)
numStimFrames = Offset - Onset;
if range(numStimFrames) > 1
    warning(['The number of stimulus frames varies by ', num2str(range(numStimFrames))])
end
numStimFrames = mean(numStimFrames);

%% Check that frames are not outside the bounds of the image
numPreFrames = settings.preTime * settings.frameRate;
numPostFrames = settings.postTime * settings.frameRate;

badEpoch = find(any(...
Onset - numPreFrames < 0 |... 
Onset + numStimFrames + numPostFrames > length(stimTrace) ...
,2));

if badEpoch
    warning(['Had to remove epoch ' mat2str(badEpoch)...
        ' because the pre or post stime stretched outside the bounds of the image.'])
    Onset(badEpoch,:) = [];
    Offset(badEpoch,:) = [];
end

%% Save important info to the output Struct
outputStruct.numStimFrames = numStimFrames;
outputStruct.stimOnsetFrames = Onset;
outputStruct.stimOffsetFrames = Offset;
outputStruct.numEpochs = length(Offset);
outputStruct.epochSeconds = settings.preTime + settings.postTime + numStimFrames/settings.frameRate; %the length of the experiment in seconds
outputStruct.frameTime = (0:length(stimTrace)-1)/settings.frameRate;

%% Plot Onset and Offset detection
figure(2)
subplot(2,1,1)
plot(Deriv)
hold on
scatter(Onset, Deriv(Onset))
scatter(Offset, Deriv(Offset))
plot([0,length(Deriv)], [threshold, threshold])
plot([0,length(Deriv)], [-1*threshold, -1*threshold])
ylabel('Derivative of pixel intensities')
hold off

subplot(2,1,2)
plot(stimTrace)
hold on
scatter(Onset, stimTrace(Onset))
scatter(Offset, stimTrace(Offset))
xlabel('zframes')
ylabel('pixelIntensity')
hold off
end