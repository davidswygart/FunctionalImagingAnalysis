function outputStruct = identifyEpochFrames(SI, settings, resampleHz)

%% Convert to 500hz remapping
remap = resampleHz/settings.frameRate;
Onset = round(SI.stimOnsetFrames * remap);
numStimFrames = round(SI.numStimFrames * remap);
numPreFrames = round(resampleHz*settings.preTime);
numPostFrames = round(resampleHz*settings.postTime);
traceLength = round(remap*length(SI.frameTime));

%% Find preFrame, Stimulus, and PostFrame Indices for 500hz
preFrames = nan(SI.numEpochs, numPreFrames);
stimFrames = nan(SI.numEpochs, numStimFrames);
postFrames = nan(SI.numEpochs, numPostFrames);
for e = 1:SI.numEpochs
    preFrames(e,:) = Onset(e)-numPreFrames:Onset(e)-1;
    stimFrames(e,:) = Onset(e):Onset(e)+numStimFrames-1;
    postFrames(e,:) = Onset(e)+numStimFrames:Onset(e)+numStimFrames+numPostFrames-1;
end
allFrames = cat(2,preFrames, stimFrames, postFrames);

%% Create times stamps for eopoch frames
numEpochFrames = length(allFrames(1,:));
epochTime = length(allFrames(1,:))/resampleHz; %Length of time in seconds for each epoch
epochFrameTime = linspace(0, epochTime, length(allFrames(1,:)));



%% Save important data to the output structure
outputStruct.traceLength = traceLength;
outputStruct.frameTime = linspace(0,SI.frameTime(end),traceLength);

outputStruct.epochTime = epochTime;
outputStruct.epochFrameTime = epochFrameTime;

outputStruct.numEpochFrames = length(allFrames(1,:));
outputStruct.allFrames = allFrames;

outputStruct.numPreFrames = length(preFrames(1,:));
outputStruct.preFrames = preFrames;

outputStruct.numStimFrames = length(stimFrames(1,:));
outputStruct.stimFrames = stimFrames;

outputStruct.numPostFrames = length(postFrames(1,:));
outputStruct.postFrames = postFrames;


end