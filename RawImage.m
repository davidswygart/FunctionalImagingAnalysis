classdef RawImage < handle
    properties
        img
        SI
        time
    end
    methods
        function obj = RawImage(rawImage_path)
            [timg, res, md, ts,pos] = fastLoadTiff(rawImage_path); 
            obj.img = permute(timg, [2, 1, 4, 3]);
            mdStruct = parseScanImageMetaData(md);
            obj.SI = mdStruct.SI;   
        end
        function createTimeImage(obj)
             fR = obj.SI.hRoiManager.scanFrameRate;
            [nR, nC, nF, ~] = size(obj.img);
            obj.time = createTimeImage(nR, nC, nF, fR);
        end
    end 
end