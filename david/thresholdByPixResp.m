%% display response by pixel
clf
resp = spot.pix.dPrime;
%resp = spot.pix.dFoF;

imshow3D(resp)
colorbar
axis image;

r = squeeze(mean(resp,[1,2], 'omitnan'))
%% threshold off of best spot size from above and overwrite pervious threshold
[~,i] = max(r);
display(['using spot size ' num2str(round(epochs.spotSizes(i)))])

bestSpot = resp(:,:,i);
sortedVals = sort(bestSpot(:));
threshVal = sortedVals(round(length(sortedVals)*.95))

raw.meetsThresh = bestSpot > threshVal;

%% Plot the thresholded pixels
imagesc(bestSpot, 'AlphaData', raw.meetsThresh)
%imagesc(bestSpot)
colorbar
axis image;