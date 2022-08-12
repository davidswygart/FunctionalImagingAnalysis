classdef ROI    
    properties
        version
        type
        count
        bbox
        coords
    end
    
    methods
        function obj = ROI(fname)
            f = fopen(fname,'rb','b'); %big endian byte order
            h = fread(f, 4, 'char');
            if ~all(char(h') == 'Iout')
                error('Not an ImageJ ROI file');
            end
            obj.version = fread(f, 1, 'uint16');
            type = fread(f,2,'uint8');
            obj.type = ROIType(type(1));
            
            obj.bbox = fread(f,4,'uint16') + 1; %one-based
            obj.count = fread(f,1,'uint16');

            if obj.type == ROIType.Line || obj.type == ROIType.Rectangle || obj.type == ROIType.Oval
                error('ROIs not yet implemented for type %s', obj.type);
            end

            fseek(f,64,'bof'); %skip the rest of the header, for now

            obj.coords = obj.bbox([2,1])' + fread(f,[obj.count,2],'ushort');
            fclose(f);
        end

        function res = apply(self, fn, img, xq, yq)
            if nargin<5
                [xq,yq] = meshgrid(1:size(img,1), 1:size(img,2));
            end
            in = inpolygon(xq, yq, self.coords(:,1),self.coords(:,2));
            
            img_ = reshape(img,numel(in),[]);
            img_ = img_(in(:),:);
            res = fn(img_);
        end        
    end
end

