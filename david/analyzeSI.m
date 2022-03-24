function [siMed, medShufSI] = analyzeSI(small, big, threshImg, threshVal)
[~,allSI] = calcAndPlotSi(small,big);


%% reshape threshimage into 1D array (rows = pixels)
threshArray = reshape(threshImg,[],1);

ne = size(allSI,3);
siArray = reshape(allSI,[],ne);

%% select pixels above threshold
tInd = find(threshArray > threshVal);

siArray = siArray(tInd,:);
siMed = median(siArray,2);
%% shuffle si values accross epochs
numPerm = 500;
siShuf = repmat(siArray,1,1,numPerm);
siShuf = Shuffle(siShuf,1);
medShufSI = squeeze(median(siShuf,2));

end