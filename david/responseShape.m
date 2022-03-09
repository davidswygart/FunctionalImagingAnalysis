loadRaw
eImg = splitRawImageIntoEpochs(raw);

%% make threshold using only the good times
pixFrac = .4;
isGood = raw.time > minMaxTime(1) & raw.time < minMaxTime(2);
goodThresh = raw.thresh;
goodThresh(~isGood) = nan;
[meetsThresh, tProj] = thresholdFromTProjection(goodThresh, pixFrac);

%% get spot sizes for each epoch and mark the ones that are adapting/bleaching with -1
isGoodEpoch = squeeze(eImg.start) > minMaxTime(1) & squeeze(eImg.start) < minMaxTime(2);
if strcmp(dataSetName, 'none')
    dataSetName = length(isGoodEpoch); %If no dataset, just say the number of epochs and it will generate spot sizes
end
spotSizes = getSplitParam(cellName,dataSetName,'curSpotSize');
spotSizes(~isGoodEpoch) = -1;

%% split epochs by spot size
[greenSpots, uSpots, einds] = splitByParameter(eImg.green,spotSizes);
[timeSpots, ~, ~] = splitByParameter(eImg.time,spotSizes);

%% pull values for each spot (restricted to threshold), sort, and decimate for easier processing
allTime = cell(1,length(greenSpots));
allVals = allTime;
for i = 1:length(greenSpots)
    t = timeSpots{i};
    [~,~,nt,ne] = size(t);
    repThresh = repmat(meetsThresh,1,1,nt,ne);
    allTime{i} = t(repThresh);
    g = greenSpots{i};
    allVals{i} = g(repThresh); 
end
%% bin and average data
 edges = -1:0.01:2;
for i = 1:length(greenSpots)
    t = allTime{i};
    g = allVals{i};
    
    [~,~,loc]= histcounts(t,edges);
    
    inBounds = loc > 0;
    
    loc = loc(inBounds);
    g = g(inBounds);
    
    allVals{i} = accumarray(loc, g) ./ accumarray(loc,1);
    allTime{i} = 0.5*(edges(1:end-1)+edges(2:end));
end

%% plot
clf
hold on
for i = 1:length(allTime)   
    x = allTime{i};
    y = allVals{i};
    
    x = x(3:end-2); % the first and last points tend to be noisy
    y = y(3:end-2);
    plot(x,y)
end


plot([0 0], ylim, '--k')
plot([1 1], ylim, '--k')
legend(num2str(uSpots))

saveas(gcf,[image,'_responseProfile.png'])

%% determine the best spot
dF = nan(length(allTime),1);
for i = 1:length(allTime)
    isPre = allTime{i} < 0;
    isStim = allTime{i} > .1 & allTime{i} < 1.1;
    
    y = allVals{i};
    
    dF(i) = mean(y(isStim)) - mean(y(isPre));
end
[~,bestInd] = max(dF)
if bestInd == length(dF)
    bestInd = bestInd-1;
end
