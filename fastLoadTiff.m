function [image,res,metadata,zs,ts,pos] = fastLoadTiff(file_name_and_path)
% [image,res,metadata,ts,zs,pos] = FASTLOADTIFF(file_name_and_path)
% IMAGE : the raw image [nx -by- ny -by- nc -by- nr -by- nz -by- nt]
% RES : (x,y,z) resolution in microns
% METADATA : the full SI metadata string, which can be converted into a
%   struct using eval()
% ZS : the actual z position of each frame, in microns [nr -by- nz -by- nt]
% TS : the actual time stamp of each frame, in seconds [nr -by- nz -by- nt]
% POS : starting (x,y,z) position of the stage, in microns
% NOTE: only works on raw scanimage files

f = fopen(file_name_and_path,'r');
header = fread(f,16,'uint8=>uint8');

if ~all(header(1:2) == uint8([73;73]))
    fclose(f);
    error('File is not saved in little endian order!!');
    %all the bytes are flipped around and need to be reversed!
end
if ~all(header(3:4) == uint8([43;0])) || typecast(header(5:8), 'uint32') ~=8
    fclose(f);
    error('File is not a valid bigTIFF!');
end
nread = typecast(header(9:16),'uint64');
skipped = fread(f,nread - 16,'uint8=>uint8'); 

%read first ifd
ntags = fread(f,1,'uint64=>uint64');
tags = fread(f,[20,ntags], 'uint8=>uint8');
next = fread(f,1,'uint64=>uint64');
ifd_size = 20*ntags + 16;
nread = nread + ifd_size;


%TODO: maybe don't hard code these?
width = typecast(tags(13:14,1),'uint16');
height = typecast(tags(13:14,2),'uint16');
bitdepth = typecast(tags(13:14,3),'uint16');
bitdepth_str = sprintf('int%d', bitdepth);
bitdepth_load = sprintf('%s=>%s', bitdepth_str, bitdepth_str);
bytesPerPixel = double(bitdepth) / 8;
bytesPerStrip = bytesPerPixel * double(width) * double(height);
loc = typecast(tags(13:end,7),'uint64');

res = double(typecast(reshape(tags(13:end,12:13),[],1),'uint32'));
res = [res(2)/res(1) res(4)/res(3)] * 1e4; %in microns


loc_loc = (8 + 20*6 + 13):(8 + 20*6 + 13 + 7); %should be 141:148
next_loc = ifd_size-7 : ifd_size;

%col 6 -> "imagedescription", col 16 -> "software" col 17 -> "artist"
% [info, nread] = read_tag_from_pointer(f, tags(:,6), nread); %frame timing
% [info, nread] = read_tag_from_pointer(f, tags(:,17), nread); %rois
[metadata, nread] = read_tag_from_pointer(f, tags(:,16), nread); %SI state
metadata = metadata';
[~,tok,~] = regexp(metadata, 'SI.hChannels.channelSave = ([^\n\r]+)','match','tokens','tokenExtents');
nchans = length(str2num(tok{1}{1})); %#ok<ST2NM>

[~,tok,~] = regexp(metadata, 'SI.hMotors.samplePosition = ([^\n\r]+)','match','tokens','tokenExtents');
pos = str2num(tok{1}{1}); %#ok<ST2NM>
%TODO: check the starting galvo position, and adjust pos so that it always
%refers to the same pixel (either center or top left of the image)

[~,tok,~] = regexp(metadata, 'SI.hStackManager.actualNumSlices = ([^\n\r]+)','match','tokens','tokenExtents');
nSlices = str2double(tok{1}{1});

%TODO: whats the deal with "zsAllActuators"?
[~,tok,~] = regexp(metadata, 'SI.hStackManager.zs = ([^\n\r]+)','match','tokens','tokenExtents');
zs = reshape(str2num(tok{1}{1}),[],1); %#ok<ST2NM>

if nSlices > 1
% [~,tok,~] = regexp(metadata, 'SI.hStackManager.actualStackZStepSize = ([^\n\r]+)','match','tokens','tokenExtents');
% res(3) = abs(str2double(tok{1}{1}));

[~,tok,~] = regexp(metadata, 'SI.hStackManager.stackZStartPos = ([^\n\r]+)','match','tokens','tokenExtents');
start = reshape(str2num(tok{1}{1}),[],1); %#ok<ST2NM>

[~,tok,~] = regexp(metadata, 'SI.hStackManager.stackZEndPos = ([^\n\r]+)','match','tokens','tokenExtents');
fin = reshape(str2num(tok{1}{1}),[],1); %#ok<ST2NM>

res(3) = abs(fin - start) / nSlices;
else 
    res(3) = nan;
end

[~,tok,~] = regexp(metadata, 'SI.hStackManager.framesPerSlice = ([^\n\r]+)','match','tokens','tokenExtents');
nReps = str2double(tok{1}{1});
nVols = [];
if nReps == Inf
    nReps = [];
    nVols = 1;
    nPer = nchans * nSlices * nVols;
else
    nPer = nchans * nSlices * nReps;
end

%TODO: correct for case where image is taken upside down? Would probably
%get the sign of the step size and then flip image along the z axis...


bytesPerFrame = bytesPerStrip + ifd_size;

file_struct = dir(file_name_and_path);
file_bytes = file_struct.bytes;
maxFrames = floor((file_bytes-nread)/ bytesPerFrame);
image = zeros(width,height,maxFrames,bitdepth_str); %TODO: uninit

if nargout>3
    ts = zeros(1,maxFrames);
    %get epoch from imagedesc
    desc_loc = (109:128)';
    [desc, nread] = read_tag_from_pointer(f, tags(:,6), nread);
    [~,tok,~] = regexp(desc', 'epoch = (\[([\d\,\.]+)\])','match','tokens','tokenExtents');
    epoch = str2num(tok{1}{1});
    
    for i = 1:maxFrames
        fread(f, loc-nread, 'uint8=>uint8'); %TODO: probably faster to use fseek?
        image(:,:,i) = fread(f, [width, height], bitdepth_load);

        if next > file_bytes || next < nread
            break
        end

        fread(f, next - loc - bytesPerStrip, 'uint8=>uint8');
        ifd = fread(f, ifd_size, 'uint8=>uint8');
%         tags = reshape(ifd(9:end-8),20,[]);
        nread = next + ifd_size;
        [desc, nread] = read_tag_from_pointer(f, ifd(desc_loc), nread);
        [~,tok,~] = regexp(desc', 'frameTimestamps_sec = (\d+\.\d+)','match','tokens','tokenExtents');
        ts(i) = str2num(tok{1}{1});
    
        loc = typecast(ifd(loc_loc),'uint64');
        next = typecast(ifd(next_loc), 'uint64');
    end
else

    for i = 1:maxFrames
        fread(f, loc-nread, 'uint8=>uint8'); %TODO: probably faster to use fseek?
        image(:,:,i) = fread(f, [width, height], bitdepth_load);

        if next > file_bytes || next < nread
            break
        end

        fread(f, next - loc - bytesPerStrip, 'uint8=>uint8');
        ifd = fread(f, ifd_size, 'uint8=>uint8');
%         tags = reshape(ifd(9:end-8),20,[]);
        nread = next + ifd_size;
%         [desc, nread] = read_tag_from_pointer(f, tags(:,6), nread);
        loc = typecast(ifd(loc_loc),'uint64');
        next = typecast(ifd(next_loc), 'uint64');
    end
end

if mod(i, nPer) ~= 0
   if nSlices == 1 && isempty(nVols) %this is okay
        nReps = i/nchans; % floor division, because i is uint64
   end
end
%TODO: will error if imaging did not run to completion except above case

image = reshape(image(:,:,1:i), width, height, nchans, nReps, nSlices, nVols);
if nargout > 3
    ts = reshape(posixtime(datetime(epoch) + seconds(ts(1:nchans:i))), nReps, nSlices, nVols);
    if nargout > 4 
        if numel(zs) == 1  %only 1 slice
            zs = repmat(zs,size(ts));
        else 
            zs = reshape(zs, nReps, nSlices, nVols);
        end
    end
end
    
end

function [tag_data, nread] =read_tag_from_pointer(f, tag, nread)

    count = typecast(tag(5:12),'uint64');
    loc = typecast(tag(13:end),'uint64');
%     fread(f, loc - nread,'uint8=>uint8');
    fseek(f, loc, -1);
    
    rat = false;
    
    switch typecast(tag(3:4),'uint16')
        case 1
            read_str = 'uint8=>uint8';
        case 2
            read_str = 'char=>char';
        case 3
            read_str = 'uint16=>uint16';
        case 5
            read_str = 'uint64=>uint64';
            rat = true;
        case 16
            read_str = 'uint64=>uint64';
    end
    
    tag_data = fread(f, count, read_str);
    
    if rat
        error('TODO: need to define behavior for rationals');
    end

    nread = loc + count;
end

% NOTE: the tag data is formatted as below
% ids = typecast(reshape(tags(1:2,:),[],1),'uint16');
% types = typecast(reshape(tags(3:4,:),[],1),'uint16');
% count = typecast(reshape(tags(5:12,:),[],1),'uint64');
% locs = typecast(reshape(tags(13:end,:),[],1),'uint64');
%
% types defines the datatype of each value
%
% if count*sizeof[types] is <8, the tag data is inlined into the next 8 bytes
% else the location is given in the next 8 bytes
%
% the tag ID defines the property identity. some are standardized, like
% image width and height. The IDs are in ascending order, but they really
% should be checked. We are ignoring that step here, assuming that they are
% well behaved across scanimage files.
%
% furthermore each ifd could contain any number of tags, but we assume
% they're the same
% 
% locs(1:2) <- the width/height of the image in pixels
% locs(3) <- the bit depth
% locs(7) <- the byte offset, from the beginning of the file, to the image
% locs(12:13) <- the byte offset of the xy-resolution, stored as a ratio



%example ImageDescription tag:
% frameNumbers = 1
% acquisitionNumbers = 1
% frameNumberAcquisition = 1
% frameTimestamps_sec = 0.000000
% acqTriggerTimestamps_sec = 0.000000
% nextFileMarkerTimestamps_sec = 0.000000
% endOfAcquisition = 0
% endOfAcquisitionMode = 0
% dcOverVoltage = 0
% epoch = [2021,6,3,15,6,47.401]
