function [paramVals] = getSplitParam(cellName, dataSetName, splitParam)
%% load cell data and get values for the split parameter
CELL_DATA_FOLDER = getenv('CELL_DATA_FOLDER');

if exist([CELL_DATA_FOLDER, filesep, cellName],'file')
    load([CELL_DATA_FOLDER, filesep, cellName],'cellData');
    
    %cellData.savedDataSets.keys;
    datasetEpochs = cellData.savedDataSets(dataSetName)';
    paramVals = cellData.getEpochVals(splitParam,datasetEpochs)';

else
    warning('no cell data found, assuming 2 spot sizes (200/1200)')
    paramVals = nan(dataSetName,1);
    paramVals(1:2:end) = 200;
    paramVals(2:2:end) = 1200;
end
end




