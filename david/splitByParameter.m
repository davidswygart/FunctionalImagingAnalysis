function [imgArray, uniqueParamVals, epochInds] = splitByParameter(eImg, paramVals)
%% assumes a struct with 4D epoch image (row,column,time,epoch number)
%% imgArray: a cell array of the image split by the parameter
%% paramVals: 1D array of each value of the split parameter
%% eInds:  1D array indicating which cell each epoch was sorted into
[nr, nc, nt, ne] = size(eImg);

if length(paramVals) > ne
    warning('there are more cellDataset epochs than image epochs')
    paramVals = paramVals(1:ne); %cut off the excess epoch parameter values
elseif length(paramVals) < ne
    warning('there are more image epochs than cellDataset epochs')
end

%% list the epoch indices corresponding to each size (in cell array)
uniqueParamVals = unique(paramVals);
imgArray = cell(1, length(uniqueParamVals));
epochInds = nan(1, length(paramVals));


for i = 1:length(uniqueParamVals)
    match = paramVals == uniqueParamVals(i);
    imgArray{i} = squeeze(eImg(:,:,:,match));
    epochInds(match) = i;
end
end




