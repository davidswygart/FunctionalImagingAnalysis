[img,res,md,zs,ts,pos] = fastLoadTiff('E:\psf_3d\region2_00005.tif');
[nx,ny,~,~,nz] = size(img);

%take median across repeats
med = double(squeeze(median(img(:,:,1,:,:), 4)));


% find peaks by gradient descent
radius = .2; %radius in microns

[Gx,Gy,Gz] = imgradientxyz(imgaussfilt3(med, radius./res));
map = zeros(size(med)); %each pixel will map to its local prominence
prom = zeros(size(med)); %we will store how many pixels map to each pixel


i = find(~map(:),1);
while ~isempty(i)
    [x,y,z] = ind2sub([nx,ny,nz],i);
    
    pix = i;
        
    while Gx(i) >= .5 || Gy(i) >= .5 || Gz(i) >= .5
        %step to the new location
        fx = max(min(round(x + sign(Gx(i))), nx), 1);
        y = max(min(round(y + sign(Gy(i))), ny), 1);
        z = max(min(round(z + sign(Gz(i))), nz), 1);
        ni = sub2ind([nx,ny,nz],x,y,z);    
        if i == ni || ismember(i,pix)
            break
        end
        i = ni;
        
        %add this location to our list...
        pix = cat(1,pix,i);
    end
    
    %pixels in our list will map to this location
    map(pix) = i;
    %store the number of mapping pixels
    prom(i) = numel(pix);    
    
    i = find(~map(:),1);
end
