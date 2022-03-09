loadRaw
eImg = splitRawImageIntoEpochs(raw);
%% get rid of epochs that occur before adaptation or after bleaching
goodE = squeeze(eImg.start) > minMaxTime(1) & squeeze(eImg.start) < minMaxTime(2);
fNames = fieldnames(eImg);
for i=1:length(fNames)
temp = eImg.(fNames{i});
eImg.(fNames{i}) = temp(:,:,:,goodE);
end

%% split by spot size
[greenSpots, spotSizes, einds] = splitByParameter(eImg.green,cellName,dataSetName,'curSpotSize');
[timeSpots, ~, ~] = splitByParameter(eImg.time,cellName,dataSetName,'curSpotSize');
[eStart, ~, ~] = splitByParameter(eImg.start,cellName,dataSetName,'curSpotSize');

%% make threshold
figure(1)
pixFrac = .4;
[meetsThresh, tProj] = thresholdFromTProjection(raw.thresh, pixFrac);
saveas(gcf,[image,'_Threshold.png'])
%% calculate the response for each epoch pixel and then average accross space
% response for each pixel
resp = cell(1,length(greenSpots));
for i = 1:length(greenSpots)
    resp{i} = epochPixelResponse(greenSpots{i}, timeSpots{i});
end

% average over space
for i = 1:length(resp)
    r = resp{i};
    temp = r.dFoF;
    repThresh = repmat(meetsThresh,1,1,size(temp,3));
    temp(~repThresh) = nan;
    r.dFoF = squeeze(mean(temp, [1,2], 'omitnan'));
    
    temp = r.dPrime;
    temp(~repThresh) = nan;
    r.dPrime = squeeze(mean(temp, [1,2], 'omitnan')); 
    resp{i} = r;
end


% plot Dprime
figure(2)
clf
subplot(2,1,1)
hold on
for i = 1:length(resp)  
    t = eStart{i};
    r = resp{i};
    y = r.dPrime;
    plot(t,y)
end
legend(num2str(round(spotSizes)))
title('pixels averaged')
xlabel('time (s)')
ylabel('dPrime')
% plot DFoF
subplot(2,1,2)
hold on
for i = 1:length(resp)  
    t = eStart{i};
    r = resp{i};
    y = r.dFoF;
    plot(t,y)
end
legend(num2str(round(spotSizes)))
xlabel('time (s)')
ylabel('dFoF')
saveas(gcf,[image,'_AdaptationBleaching.png'])

%% calculate the response for each for each frame
% 
% 
% %% calculate values for each spot
% spotsR = struct;
% for i = 1:length(eImg.uSpots)
%     inds = eImg.spotInd{i};
%     
%     for f = 1:length(fNames)
%         
%     end
%     
%     % copied from epochs
%     spotsR(i).pre = eImg.pre(:,:,inds);
%     spotsR(i).stim = eImg.stim(:,:,inds);
%     spotsR(i).dF = eImg.dF(:,:,inds);
%     spotsR(i).dFoF = eImg.dFoF(:,:,inds);
%     spotsR(i).dPrime = eImg.dPrime(:,:,inds); 
%     
%     % calculate means
%     spotsR(i).avgPre = mean(spotsR(i).pre,3);
%     spotsR(i).avgStim = mean(spotsR(i).stim,3);
%     spotsR(i).avgDF = mean(spotsR(i).dF,3);
%     spotsR(i).avgDFoF = mean(spotsR(i).dFoF,3);
%     spotsR(i).avgDPrime = mean(spotsR(i).dPrime,3);
%     
%     % calculate means
%     spotsR(i).stdPre = std(spotsR(i).pre, 0,3);
%     spotsR(i).stdStim = std(spotsR(i).stim, 0,3);
%     spotsR(i).stdDF = std(spotsR(i).dF,0 ,3);
%     spotsR(i).stdDFoF = std(spotsR(i).dFoF, 0,3);
%     spotsR(i).stdDPrime = std(spotsR(i).dPrime, 0,3);
% end
% 
% %% plot SMS (averaged accross epochs then averaged accross time)
% clf
% x = eImg.uSpots;
% y = nan(1,length(spotsR));
% err = nan(1,length(spotsR));
% for i = 1:length(spotsR)
%     temp = spotsR(i).avgDFoF(raw.meetsThresh);
%     %temp = spot(i).avgDPrime(raw.meetsThresh);
%     %temp = spot(i).avgDF(raw.meetsThresh);
%     %temp = spot(i).avgPre(raw.meetsThresh);
%     %temp = spot(i).avgStim(raw.meetsThresh);
%     
%     y(i) = mean(temp);
%     err(i) = std(temp);
% end
% 
% % repThresh = repmat(raw.meetsThresh,1,1,length(x));
% % y = mean(spot.pix.dFoF(repThresh));
% %y = spot.fram.dFoF;
% %y = spot.fram.avgPre;
% %y = spot.fram.avgStim;
% %y = spot.fram.stdPre;
% %y = spot.fram.stdStim;
% 
% errorbar(x,y,err)
% xlabel('spot diameter (um)')
% ylabel('dPrime')
% 
% hold on
% plot([0,1300], [0 0], '--k') % plot y = 0
% %plot([0,1300], [2 2], '--k') % plot y = 2 standard deviations
% 
% 
% [~,bestInd] = max(y);
% %% plot signal vs. error scatter for each pixel that meets threshold
% x = spotsR(bestInd).stdDPrime;
% y = spotsR(bestInd).avgDPrime;
% % x = spot(bestInd).stdDFoF;
% % y = spot(bestInd).avgDFoF;
% %  x = spot(bestInd).stdDF;
% %  y = spot(bestInd).avgDF;
% 
% % x = x(:);
% % y = y(:);
% 
% x = x(raw.meetsThresh);
% y = y(raw.meetsThresh);
% 
% clf
% hold on
% scatter(x,y)
% 
% plot([0 max([x;y])],[0 max([x;y])]) %unity line
% 
% ylabel('signal')
% xlabel('error')
% %% plot each pixel dFoF over time
% % clf
% % hold on
% % %for i = 1:spot.n
% % i = bestSpotInd;
% % inds = epochs.spotInd{i};
% % 
% % t = epochs.startTime(inds);
% % img = epochs.pix.dFoF(:,:,inds);
% % 
% % [ri ci] = find(raw.meetsThresh);
% % 
% % for ii = 1:length(ri)
% %     y = squeeze(img(ri(ii),ci(ii),:));
% %     plot(t,y)
% % end
% % %end
% % 
% % %scatter(minMaxTime, [0,0], 'filled','k')
% % 
% % %legend(num2str(round(epochs.uSpots)))
% % title('response over time')
% % xlabel('time (s)')
% % ylabel('response')
% % 
% % 
