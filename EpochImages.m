classdef EpochImages
    properties
        img5D
        dimensions = ['lines (Y)','columns (X)','time (frames)','Epochs', 'Channel'];
        channelOrder
        cellName
        dataSetName
        pxOn
        rate
        additionalMetaData
        figMap = containers.Map({'plotRaw','splitByParam', 'appliedMask','overTrials'}, 1:4);
        
    end
    
    methods (Hidden = true) %% Building functions
        function obj = EpochImages(rawImage, cellName, dataSetName)       
           % Set image properties
           obj.cellName = cellName;
           obj.dataSetName = dataSetName;
           obj.pxOn = rawImage.pxOn;
           obj.rate = rawImage.rate;
           obj.additionalMetaData = rawImage.additionalMetaData;
           obj.verifyMatchingEpochs();
           

           
           % Build Epoch Image
            chOrd = rawImage.channelOrder;
            chOrd(chOrd == "pixOrder") = "time";
            obj.channelOrder = chOrd;
            obj.img5D = obj.buildEpochImg(rawImage.img);
            
        end
        function verifyMatchingEpochs(obj)
           %Check that Symphony and the image have same number of epochs
           numSymphonyEpochs = length(obj.epochsByParam('preTime'));
           numImageEpochs = length(obj.pxOn);
           if (numImageEpochs ~= numSymphonyEpochs)
               msg = 'the number of image triggers and dataset epochs are not all equal';
               msg2 = sprintf('image trigs = %d', numImageEpochs);
               msg3 = sprintf('dataset epochs = %d', numSymphonyEpochs);
               error([msg newline msg2 newline msg3]);
           end
        end
        
        function epochImg = buildEpochImg(obj, img)
            frameRate = obj.rate.frame;
            [~, preTime] = obj.epochsByParam('preTime');
            [~, stimTime] = obj.epochsByParam('stimTime');
            [~, tailTime] = obj.epochsByParam('tailTime');
            
            if length([preTime,stimTime,tailTime]) > 3
                error('more than 1 length of preTime/stimTime/tailTime detected for epochs')
            end
            
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
            parse(p,varargin{:});
            
            if isempty(p.Results.ChannelNames)
                ch = 1:sz(5);
            else
                [~,ch] = ismember(p.Results.ChannelNames,obj.channelOrder);
            end

            subImage = obj.img5D(:,:,p.Results.Frames,p.Results.Epochs,ch);
            subImage = obj.apply2DMask(subImage,p.Results.xyMask);
        end
        
        function [SNR, dFoF, avgStim, avgPre] = epochResponses(obj,varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'channel',"green");
            addParameter(p,'ratioChannel',[]);
            addParameter(p,'method','SNR'); % 'SNR', 'dF/F'
            parse(p,varargin{:});
            
            
            
            chs = ["time",p.Results.channel,p.Results.ratioChannel]; %all channel names
            img = obj.getSubImage(p.Unmatched,'ChannelNames',chs);
            
            t = img(:,:,:,:,1);
            if isempty(p.Results.ratioChannel)
                v = img(:,:,:,:,2);
            else
                v = img(:,:,:,:,2) ./ img(:,:,:,:,3); %divide by 3rd channel (red)
            end
            
            preVals = v; %copy all values from analysis channel
            preVals(t > 0) = nan; %nan out values that are not pre-Stim
            avgPre = squeeze(mean(preVals,3,'omitnan')); %average accross frames (time) and squeeze out the time dimension
            
            stimVals = v; %copy all values from analysis channel
            [~, stimTime] = obj.epochsByParam('stimTime');
            stimInds = (t > 0) & (t < (stimTime / 1000));           
            stimVals(~stimInds) = nan; %nan out values that are not StimTime
            avgStim = squeeze(mean(stimVals,3,'omitnan')); %average accross frames (time) and squeeze out the time dimension
            
            stdStim = squeeze(std(stimVals,0,3,'omitnan'));           
            dF = avgStim - avgPre;
            
            %if strcmp(p.Results.method,'SNR')
                SNR = dF ./ stdStim;
            %else if strcmp(p.Results.method,'dF/F')
                dFoF = dF ./ avgPre;
            %end               
        end
        
        function image = apply2DMask(obj,image,mask)
            if ~isempty(mask)
                sz = size(image);
                expandedMask = repmat(mask,[1,1,sz(3:end)]);
                image(~expandedMask) = nan;
                
                figure(obj.figMap('appliedMask'));
                imagesc(image(:,:,1));
                title('mask')
            end
        end
        
        function [l, errFill, fig] = plotRaw(obj,varargin)
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'channelY',"green");
            addParameter(p,'channelX',"time");
            addParameter(p,'binTime',0.05);
            addParameter(p,'figName','plotRaw');
            parse(p,varargin{:});
            
            xImg = obj.getSubImage(p.Unmatched,'ChannelNames',p.Results.channelX);
            yImg = obj.getSubImage(p.Unmatched,'ChannelNames',p.Results.channelY);

%           binscatter(xImg(:),yImg(:), 240)


            x_all = xImg(~isnan(xImg));
            y_all = yImg(~isnan(yImg));
            edges = min(x_all):p.Results.binTime:max(x_all)+p.Results.binTime;
            binNum = discretize(x_all,edges); %bin datapoints by X and get their bin number
            avgY = accumarray(binNum,y_all,[],@(v)mean(v));
            semY = accumarray(binNum,y_all,[],@(v)std(v)./sqrt(length(v)));
           
            xmid = 0.5*(edges(1:end-1)+edges(2:end))';
            
            
            fig = figure(obj.figMap(p.Results.figName));
            if strcmp(p.Results.figName,'plotRaw')
                fig.Name = 'plotRaw';
                clf
                hold on
            end
            
            l = plot(xmid,avgY);
            errFill = fill([xmid; flip(xmid)],[avgY-semY; flip(avgY+semY)],l.Color,'linestyle','none','FaceAlpha',.5);
            xlabel(p.Results.channelX);
            ylabel(p.Results.channelY);
        end
        
        function fig = plotSplitByParam(obj,param,varargin)
            [eInd, uniqueVals] = obj.epochsByParam(param);
            
            fig = figure(obj.figMap('splitByParam'));
            fig.Name = 'splitByParam';
            clf
            hold on
            lines = matlab.graphics.chart.primitive.Line.empty(0,length(uniqueVals));
            for i = 1:length(uniqueVals)
                [lines(i),~,~] = obj.plotRaw(varargin{:}, 'Epochs',i == eInd, 'figName', 'splitByParam'); 
            end
            title(sprintf('split by %s', param));
            legend(lines,num2str(uniqueVals'));
        end
        
        function fig = plotOverTrials(obj,varargin) 
            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p,'paramNameAndValue',[]); %requires both a name and value in cell
            parse(p,varargin{:});
            
            if length(p.Results.paramNameAndValue) == 2
                [e, ~] = obj.epochsByParam(p.Results.paramNameAndValue{:});
                eInds = find(e);
            else
                numEpochs = size(obj.img5D, 4);
                eInds = 1:numEpochs;
            end
                       
            fig = figure(obj.figMap('overTrials'));
            fig.Name = 'overTrials';
            clf
            hold on
            lines = matlab.graphics.chart.primitive.Line.empty(0,length(eInds));
            for i = 1:length(eInds)
                [lines(i),~,~] = obj.plotRaw(p.Unmatched, 'Epochs',eInds(i), 'figName', 'overTrials'); 
            end
            title('Response by epoch');
            legend(lines,num2str(eInds'));        
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