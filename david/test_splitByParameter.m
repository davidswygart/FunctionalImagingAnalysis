function [imgArray, uniqueParamVals, epochInds] = splitByParameter(eImg, cellName, dataSetName, splitParam)
%% assumes a struct with 4D epoch image (row,column,time,epoch number)
%% imgArray: a cell array of the image split by the parameter
%% paramVals: 1D array of each value of the split parameter
%% eInds:  1D array indicating which cell each epoch was sorted into

[nr, nc, nt, ne] = size(eImg);

%% load cell data and get values for the split parameter
CELL_DATA_FOLDER = getenv('CELL_DATA_FOLDER');

if exist([CELL_DATA_FOLDER, filesep, cellName],'file')
    load([CELL_DATA_FOLDER, filesep, cellName],'cellData');
    
    %cellData.savedDataSets.keys;
    datasetEpochs = cellData.savedDataSets(dataSetName)';
    epochParamVals = cellData.getEpochVals(splitParam,datasetEpochs)';
    %epochs.spotSizes = spotSizes(1:numImgEpochs);
    if length(epochParamVals) > ne
        warning('there are more cellDataset epochs than image epochs')
        epochParamVals = epochParamVals(1:ne); %cut off the excess epoch parameter values
    elseif length(epochParamVals) < ne
        warning('there are more image epochs than cellDataset epochs')
    end
else
    warning('no cell data found, assuming 2 spot sizes (200/1200)')
    epochParamVals = nan(ne,1);
    epochParamVals(1:2:end) = 200;
    epochParamVals(2:2:end) = 1200;
end
%% list the epoch indices corresponding to each size (in cell array)
uniqueParamVals = unique(epochParamVals);
imgArray = cell(1, length(uniqueParamVals));
epochInds = nan(1, length(epochParamVals));

for 1:
for i = 1:length(uniqueParamVals)
    match = epochParamVals == uniqueParamVals(i);
    imgArray{i} = squeeze(eImg(:,:,:,match));
    epochInds(match) = i;
end
end




