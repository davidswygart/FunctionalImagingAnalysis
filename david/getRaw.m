function [raw,metaData]= getRaw(image)
path = 'C:\Users\david\Desktop\2p_copies';

%% load the raw image and metadata
[tempImg, res, md, ~,pos] = fastLoadTiff([path filesep image]); %not using timestamps
raw = permute(tempImg, [2, 1, 4, 3]);

% parse metatdata into struct to be more readable
metaData = parseScanImageMetaData(md);
end