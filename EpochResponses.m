classdef EpochResponses
    properties
        cellName
        dataSetName
        img
        dimensions = ['lines (Y)','columns (X)','Epochs'];
        figMap = containers.Map({'blank'}, 1:1);
        
    end
    
    methods (Hidden = true) %% Building functions
        function obj = EpochResponses(img, cellName, dataSetName)       
           % Set image properties
           obj.img = img;
           obj.cellName = cellName;
           obj.dataSetName = dataSetName;
        end
    end
    methods
        function [eInd, uniqueVals] = epochsByParam(obj,paramName, varargin) %varargin is optional value for the parameter name
            CELL_DATA_FOLDER = getenv('CELL_DATA_FOLDER');
            load([CELL_DATA_FOLDER, filesep, obj.cellName],'cellData');
            
            try
                datasetEpochs = cellData.savedDataSets(obj.dataSetName);
            catch
                names = cellData.savedDataSets.keys;
                errMsg = sprintf('Unable to find Dataset named %s\n Available Datasets for %s:\n', obj.dataSetName, obj.cellName);
                nameMsg = sprintf('%s \n', names{:});
                error([errMsg,nameMsg]);
            end
            
            paramVals = cellData.getEpochVals(paramName,datasetEpochs);
            
            if isempty(varargin)
                uniqueVals = unique(paramVals);
            else
                uniqueVals = varargin{:};
            end
            
            [~,eInd] = ismembertol(paramVals,uniqueVals,.001); %use .001 as tolerance because of strangeness in comparing doubles
        end
    end
end