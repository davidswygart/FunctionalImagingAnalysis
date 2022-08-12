classdef ROIs  < handle
    properties
        image
        roiData = struct.empty;
        axTrace = nan;
        axResp = nan;
        lsResp = [];
        
        plotFunc = [];
        plotParams = {};
        
    end
    
    
    methods
        function obj = ROIs(varargin)
            obj.image = gca;

%            obj.createTraceImage();
%obj.createRespImage();
            if ~isempty(varargin)
                obj.plotFunc = varargin{1};
                if length(varargin) > 1
                    obj.plotParams = varargin{2};
                end
            end

            obj.addROI();%add the first ROI
       end
        
        function createTraceImage(obj)
            %figure(111)
            obj.axTrace = gca;
            
            %imagesc(obj.image, limit)
            %scnan(obj.image, obj.thresh)
            %caxis([0 1])
            %colorbar
            

            for i = 1:length(obj.roiData)
                roi = obj.createROI(obj.roiData(i)); %recreate the roi
                obj.roiData(i).handle = roi; %update the roi handle
            end
        end
        
        function createRespImage(obj)
            figure(112);
            obj.axResp = gca;
            
            for i = 1:length(obj.roiData)
                if ~ishandle(obj.roiData(i).handle) %verify the handle is still valid
                    roi = obj.createROI(obj.roiData(i)); %recreate the roi
                    obj.roiData(i).handle = roi; %update the roi handle
                end
                obj.plotROIResp(obj.roiData(i).handle);
            end

                
                
        end
        
        function addROI(obj)
            newROI = obj.createROI();
            
            ind = length(obj.roiData) + 1;
            newROI.Label = num2str(ind);
            colors = colororder();
            cInd = mod(ind,size(colors,1))+1; %loop through the colors (starts at 2, but that isn't a problem)
            newROI.Color = colors(cInd,:); %Redo so that line chooses color
            obj.roiData(ind).position = newROI.Position;
            obj.roiData(ind).color = newROI.Color;
            obj.roiData(ind).handle = newROI;
            
            if ~isempty(obj.plotFunc)
                obj.plotROIResp(newROI);
            end
        end
        
        function roi = createROI(obj, varargin) 
            if ~isgraphics(obj.axTrace)
                obj.createTraceImage(); %remake figure if handle is no longer valid
            end
            
            if isempty(varargin)
                fprintf('Draw an ROI \n');
                roi = drawpolygon(obj.axTrace); 
            else
                roiD = varargin{1};
                c = roiD.color;
                p = roiD.position;
                roi = images.roi.Polygon(obj.axTrace, 'Color', c,'Position', p);
            end
                     
            addlistener(roi,'ROIMoved',@obj.allevents);
            addlistener(roi,'DrawingFinished',@obj.allevents);
            addlistener(roi,'DeletingROI',@obj.allevents); 
        end
        
        function plotROIResp(obj,roi)
            if ~isgraphics(obj.axResp)
                obj.createRespImage(); %axis are invalid, just recreate the whole figure.  hopefully this doesn't become an endless loop.
            else
                mask = roi.createMask();
                [l] = obj.plotFunc(obj.plotParams{:},'xyMask', mask, 'ax', obj.axResp);
                %l.Color = roi.Color;
                %err.FaceColor = roi.Color;
                obj.lsResp = [obj.lsResp, l]; 
                %label = 1:length(obj.lsResp);
                %legend(obj.lsResp, num2str(label'))
            end
        end
        
        function allevents(r,p,e)            
           % display(evt.EventName);
%             switch(evname)
%                 case{'MovingROI'}
%                     disp(['ROI moving Previous Position: ' mat2str(evt.PreviousPosition)]);
%                     disp(['ROI moving Current Position: ' mat2str(evt.CurrentPosition)]);
%                 case{'ROIMoved'}
%                     disp(['ROI moved Previous Position: ' mat2str(evt.PreviousPosition)]);
%                     disp(['ROI moved Current Position: ' mat2str(evt.CurrentPosition)]);
%             end
        end
    end
end