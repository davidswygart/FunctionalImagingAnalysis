function image = fastLoadTiffGen(file_name_and_path)
% image = FASTLOADTIFFGEN(file_name_and_path)
% IMAGE : the raw image [nx -by- ny -by- nf]
% NOTE: For generic tiff files. Does not read metadata.

f = fopen(file_name_and_path,'r');
header = fread(f,8,'uint8=>uint8');

if ~all(header(1:2) == [0x49;0x49])
    fclose(f);
    error('File is not saved in little endian order!!');
    %all the bytes are flipped around and need to be reversed!
end
if ~all(header(3:4) == [0x2A;0x00])
    fclose(f);
    error('File is not a valid tiff!');
end
prev = typecast(header(5:8),'uint32');
fseek(f,prev,'bof');

%read first ifd
ntags = fread(f,1,'uint16=>uint32');
tags = fread(f,[12,ntags], 'uint8=>uint8');
next = fread(f,1,'uint32=>uint32');

%TODO: maybe don't hard code these?
width = typecast(tags(9:12,1),'uint32');
height = typecast(tags(9:12,2),'uint32');
bitdepth = typecast(tags(9:10,3),'uint16');
bitdepth_str = sprintf('int%d', bitdepth);
bitdepth_load = sprintf('%s=>%s', bitdepth_str, bitdepth_str);
bytesPerPixel = double(bitdepth) / 8;
bytesPerStrip = bytesPerPixel * double(width) * double(height);

% fseek(f,typecast(tags(9:end,6),'uint32'),'bof');
% desc = char(fread(f,typecast(tags(5:8,6),'uint32'),'char')');

% fseek(f,typecast(tags(9:end,11),'uint32'),'bof');
% res = fread(f,4,'uint32=>double');
% res = [res(1)/res(2), res(3)/res(4)];


file_bytes = dir(file_name_and_path).bytes;
maxFrames = floor(file_bytes/ bytesPerStrip);
image = zeros(width,height,maxFrames,bitdepth_str); %TODO: uninit

    
for i = 1:maxFrames
    fseek(f,typecast(tags(9:end,7),'uint32'),'bof');
    image(:,:,i) = fread(f, [width, height], bitdepth_load);

    if next > file_bytes || next < prev
        break
    end
    fseek(f, next, 'bof');
    prev = next;
    ntags = fread(f,1,'uint16=>uint32');
    tags = fread(f,[12,ntags], 'uint8=>uint8');
    next = fread(f,1,'uint32=>uint32');
end

image = image(:,:,1:i);

end

