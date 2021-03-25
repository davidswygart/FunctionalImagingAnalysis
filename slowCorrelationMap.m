function map = slowCorrelationMap(image, sigma, radius)
padded_image = padarray(image,[radius,radius],nan);
map = zeros(size(image,1), size(image,2));

flattened = reshape(padded_image,[],size(image,3))';
stride = 2*radius + size(image,1);
ind_mat = radius + (1:size(image,1)) + stride * (radius - 1 + (1:size(image,2))');
%these are the indices of the actual image inside the flattened padded
%image

offsets = (-radius:radius) + stride*(-radius:radius)';
%weights are standard gaussian
weights = exp(-((-radius:radius).^2/(2*sigma.^2) + (-radius:radius)'.^2/(2*sigma.^2))) / sqrt((2*pi).^2 .* sigma);

for i = 1:2*radius+1
    for j= 1:2*radius+1
        if i==radius+1 && j==radius+1
            continue
        end
        response = weights(i,j) .*arrayfun(@(x) corr(flattened(:,x), flattened(:,x+offsets(i,j))), ind_mat);
        keep = ~isnan(response);
        map(keep) = map(keep) + response(keep);
    end
end
map = map';

end