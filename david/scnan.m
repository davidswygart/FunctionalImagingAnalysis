function scnan(image, threshold)
image = squeeze(image);
threshold = squeeze(threshold);

imagesc(image, 'AlphaData', threshold)
colorbar
end