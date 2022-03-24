function [paramVals] = getSplitParam2(cellName, dataSetName, splitParam)
%% load cell data and get values for the split parameter
CELL_DATA_FOLDER = getenv('CELL_DATA_FOLDER');

load([CELL_DATA_FOLDER, filesep, cellName],'cellData');

%cellData.savedDataSets.keys;
datasetEpochs = cellData.savedDataSets(dataSetName)';
paramVals = cellData.getEpochVals(splitParam,datasetEpochs)';
end




