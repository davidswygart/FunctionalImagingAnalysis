repThresh = repmat(raw.meetsThresh ,1,1,epochs.nTotalFram, epochs.n); %tProj threshold repeated

%for each frame
temp(~repThresh) = nan; %only include values that meet our previously decided threshold
fram.avgPre = squeeze(mean(temp, [1,2,3], 'omitnan'));
fram.stdPre = squeeze(std(temp,0, [1,2,3], 'omitnan'));

%for each frame
temp(~repThresh) = nan;
fram.avgStim = squeeze(mean(temp, [1,2,3], 'omitnan'));
fram.stdStim = squeeze(std(temp,0, [1,2,3], 'omitnan'));

%% calculate frame measures of response
fram.dF = fram.avgStim - fram.avgPre;
fram.dPrime = fram.dF ./ ((fram.stdStim + fram.stdPre) / 2);
fram.dFoF = fram.dF ./ fram.avgPre;