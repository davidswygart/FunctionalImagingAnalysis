function imgStruct = getRawImage(rawImage_path, channelOrder)
[timg, res, md, ts,pos] = fastLoadTiff(rawImage_path); 
img = permute(timg, [2, 1, 4, 3]);
 
[nx,ny,nz,~] = size(img);
frameRate = 1/mean(diff(ts));

if (regexp(md, 'bidirectional = true'))
    bidirectional = true;
elseif (regexp(md, 'bidirectional = false'))
    bidirectional = false;
else
    warning('Could not determine from Meta-Data if image was Bidirecitonal or Unidirectional: Assuming Bidirectional')
    bidirectional  = true;
end

% Get Pixel order and add as an additional channel
pixOrder = makePixelNums([nx,ny,nz], bidirectional);
img = cat(4,int64(img), pixOrder);
channelOrder = [channelOrder, "pixOrder"];

% Find trigger times based on image 
stimChannel = channelOrder == 'Stim';
pxOn = getTriggerTime(pixOrder, img(:,:,:,stimChannel));


imgStruct = struct;
imgStruct.img = img;
imgStruct.pxOn = pxOn;
imgStruct.channelOrder = channelOrder;
imgStruct.resolution = struct('x', res(1), 'y', res(2)); %Need to check that this is correct.
imgStruct.rate = struct('frame',frameRate, 'line',frameRate*ny,  'pixel', frameRate*nx*ny);
imgStruct.position = pos;
imgStruct.additionalMetaData = md;
imgStruct.bidirectional = bidirectional;
end