function map = fastCorrelationMap(image, sigma)
c_image = image - mean(image,3);

map = sum(imgaussfilt(c_image, sigma).^2, 3)./ imgaussfilt(sum(c_image.^2,3), sigma);
end