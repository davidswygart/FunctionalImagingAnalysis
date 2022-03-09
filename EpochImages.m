classdef EpochImages < handle
    properties
        SI
        img5D
        cellData
    end
    
    methods (Hidden = true) %% Building functions
        function obj = EpochImages(rawImage, cellName, dataSetName, triggerChannel)
            obj.SI = rawImage.SI;
            
            CELL_DATA_FOLDER = getenv('CELL_DATA_FOLDER');
            cdStruct = load([CELL_DATA_FOLDER, filesep, cellName],'cellData');
            obj.cellData = cdStruct.cellData;
            
            onLines = obj.detectStimLines(rawImage.img(:,:,:,triggerChannel));
            
            %Check that Symphony and the image have same number of epochs
            numImageEpochs = length(onLines);
            numSymphonyEpochs = length(obj.cellData.savedDataSets(dataSetName));
            if (numImageEpochs ~= numSymphonyEpochs)
                msg = 'the number of image triggers and dataset epochs are not all equal';
                msg2 = sprintf('image trigs = %d', numImageEpochs);
                msg3 = sprintf('dataset epochs = %d', numSymphonyEpochs);
                warning([msg newline msg2 newline msg3]);
            end
            
            preTime = 1000;
            stimTime = 1000;
            tailTime = 500;
            
            numPreFrames = ceil(frameRate * preTime / 1000);
            framesPerEpoch = ceil(frameRate * (preTime + stimTime + tailTime) / 1000);
            [nx,ny,nz,nCh] = size(img);
            [~,~,zs] = ind2sub([nx,ny,nz], obj.pxOn);
            startZ = zs - numPreFrames;
            endZ = startZ + framesPerEpoch-1;
            
            numEpochs = length(obj.pxOn);
            
            pxOrdCh = obj.channelOrder == "time";
            epochImg = nan(nx,ny,framesPerEpoch,nCh,numEpochs);
            
            for i=1:numEpochs
                epochImg(:,:,:,:,i) = img(:,:,startZ(i):endZ(i),:);
                epochImg(:,:,:,end,i) = ((epochImg(:,:,:,pxOrdCh,i) - obj.pxOn(i))/obj.rate.pixel); %assume last channel is pixel count
            end
            epochImg = permute(epochImg,[1,2,3,5,4]);
            
        end

    end
    methods
        function subImage = getSubImage(obj, varargin)
            sz = size(obj.img5D);
            
            p = inputParser;
            addParameter(p,'xyMask',[]);
            addParameter(p,'Frames',1:sz(3));
            addParameter(p,'Epochs',1:sz(4));
            addParameter(p,'ChannelNames',[]);
            addParameter(p,'paramNameAndValue',[]); %requires both a name and value in cell
            parse(p,varargin{:});
            
            if length(p.Results.paramNameAndValue) == 2
                [e, ~] = obj.epochsByParam(p.Results.paramNameAndValue{:});
                e = find(e);
            else
                e = p.Results.Epochs;
            end
            
            if isempty(p.Results.ChannelNames)
                ch = 1:sz(5);
            else
                [~,ch] = ismember(p.Results.ChannelNames,obj.channelOrder);
            end
            
            subImage = obj.img5D(:,:,p.Results.Frames,e,ch);
            subImage = obj.apply2DMask(subImage,p.Results.xyMask);
        end
        function [SNR, dFoF, avgStim, avgPre, dF, preVals, stimVals] = epochResponses(obj,varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'channel',"green");
            addParameter(p,'ratioChannel',[]);
            addParameter(p,'method','SNR'); % 'SNR', 'dF/F'
            parse(p,varargin{:});
            
            chs = ["time",p.Results.channel,p.Results.ratioChannel]; %all channel names
            img = obj.getSubImage(p.Unmatched,'ChannelNames',chs);
            
            
            t = img(:,:,:,:,1);
            v = img(:,:,:,:,2);
            
            %             sz = size(v);
            %             for fnum = 1:sz(3)
            %                 for enum = 1:sz(4)
            %                     filt = ones(3);
            %                     v(:,:,fnum,enum) = filter2(filt, v(:,:,fnum,enum));
            %                 end
            %             end
            
            preVals = v; %copy all values from analysis channel
            preVals(t > 0) = nan; %nan out values that are not pre-Stim
            avgPre = squeeze(mean(preVals,3,'omitnan')); %average accross frames (time) and squeeze out the time dimension
            
            stimVals = v; %copy all values from analysis channel
            [~, stimTime] = obj.epochsByParam('stimTime');
            stimInds = (t > 0) & (t < (stimTime / 1000));
            stimVals(~stimInds) = nan; %nan out values that are not StimTime
            avgStim = squeeze(mean(stimVals,3,'omitnan')); %average accross frames (time) and squeeze out the time dimension
            
            %stdStim = squeeze(std(stimVals,0,3,'omitnan'));
            stdPre = squeeze(std(stimVals,0,3,'omitnan'));
            dF = avgStim - avgPre;
            
            if ~isempty(p.Results.ratioChannel)
                meanRed = squeeze(mean(img(:,:,:,:,3), 3)); %average red through each epoch
                dF = dF ./ meanRed;
            end
            
            %if strcmp(p.Results.method,'SNR')
            SNR = dF ./ stdPre;
            %else if strcmp(p.Results.method,'dF/F')
            dFoF = dF ./ avgPre;
            %end
            
            
            
        end
        function image = apply2DMask(obj,image,mask)
            if ~isempty(mask)
                sz = size(image);
                expandedMask = repmat(mask,[1,1,sz(3:end)]);
                image(~expandedMask) = nan;
                
                currentFig = gcf;
                figure(obj.figMap('appliedMask'));
                imagesc(image(:,:,1));
                title('mask')
                figure(currentFig);
            end
        end
        function [l, errFill] = plotRaw(obj,varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'channelY',"green");
            addParameter(p,'channelX',"time");
            addParameter(p,'binTime',0.01);
            addParameter(p,'ax',[]);
            parse(p,varargin{:});
            
            xImg = obj.getSubImage(p.Unmatched,'ChannelNames',p.Results.channelX);
            yImg = obj.getSubImage(p.Unmatched,'ChannelNames',p.Results.channelY);
            
            %           binscatter(xImg(:),yImg(:), 240)
            
            x_all = xImg(~isnan(xImg));
            y_all = yImg(~isnan(yImg));
            
            startTime = min(x_all);
            endTime = max(x_all);
            numBins = (endTime - startTime) / p.Results.binTime;
            edges = linspace(startTime, endTime, numBins);
            binInd = discretize(x_all,edges); %bin datapoints by X and get their bin number
            avgY = accumarray(binInd,y_all,[],@(v)mean(v));
            xmid = 0.5*(edges(1:end-1)+edges(2:end))';
            
            semY = accumarray(binInd,y_all,[],@(v)std(v)./sqrt(length(v)));
            errX = [xmid; flip(xmid)];
            errY = [avgY-semY; flip(avgY+semY)];
            
            ax = p.Results.ax;
            if isempty(ax) %if no axis were specified, create a new figure
                fig = figure(obj.figMap('plotRaw'));
                clf
                ax = gca;
                fig.Name = 'plotRaw';
            end
            
            
            % Do actual plotting
            hold on
            l = plot(ax,xmid,avgY,'LineWidth',2);
            errFill = fill(ax,errX,errY, l.Color,'linestyle','none','FaceAlpha',.5);
            
            
            
            xlabel(ax,p.Results.channelX);
            ylabel(ax,p.Results.channelY);
        end
        function fig = plotSplitByParam(obj,param,varargin)
            [eInd, uniqueVals] = obj.epochsByParam(param);
            
            fig = figure(obj.figMap('splitByParam'));
            clf
            fig.Name = 'splitByParam';
            ax = gca;
            hold on
            lines = matlab.graphics.chart.primitive.Line.empty(0,length(uniqueVals));
            for i = 1:length(uniqueVals)
                [lines(i),~] = obj.plotRaw(varargin{:}, 'Epochs',i == eInd, 'ax',ax);
            end
            title(sprintf('split by %s', param));
            legend(lines,num2str(uniqueVals'));
        end
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





%         function fig = plotOverTrials(obj,varargin)
% %             p = inputParser;
% %             p.KeepUnmatched = true;
% %             addParameter(p,'paramNameAndValue',[]); %requires both a name and value in cell
% %             parse(p,varargin{:});
% %
% %             if length(p.Results.paramNameAndValue) == 2
% %                 [e, ~] = obj.epochsByParam(p.Results.paramNameAndValue{:});
% %                 eInds = find(e);
% %             else
% %                 numEpochs = size(obj.img5D, 4);
% %                 eInds = 1:numEpochs;
% %             end
%
%             fig = figure(obj.figMap('overTrials'));
%             fig.Name = 'overTrials';
%             clf
%             hold on
%             lines = matlab.graphics.chart.primitive.Line.empty(0,length(eInds));
%             for i = 1:length(eInds)
%                 [lines(i),~,~] = obj.plotRaw(p.Unmatched, 'Epochs',eInds(i), 'figName', 'overTrials');
%             end
%             title('Response by epoch');
%             legend(lines,num2str(eInds'));
%         end