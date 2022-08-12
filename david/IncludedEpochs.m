%% plot average green accross time
clf
hold on

%plot the green signal
yyaxis right
y = squeeze(mean(raw.green, [1,2]));
plot(raw.framT,y)

%plot the stimulus signal
yyaxis left
y = squeeze(raw.stim(1,1,:));
plot(raw.framT,y)

%show which times I am using
scatter(minMaxTime, [max(y), max(y)],'r')

%add text labels
x = raw.framT(epochs.on(:,3));
y = zeros(size(x));
text(x,y, string(round(epochs.spotSizes)))

hold off
%% estimage the good epochs
epochs.startTime = raw.framT(epochs.on(:,3));
epochs.isGood = epochs.startTime >  minMaxTime(1) & epochs.startTime < minMaxTime(2);

%% plot tProj of good times in the threshold channel
raw.isGood = raw.time > minMaxTime(1) & raw.time < minMaxTime(2);
raw.tProj = raw.thresh;
raw.tProj(~raw.isGood) = nan;
raw.tProj = mean(raw.tProj,3,'omitnan');

clf
imagesc(raw.tProj)
colorbar
axis image;
title('time projection')

%% threshold t-projection
s = sort(raw.tProj(:));
thresh = s(round(length(s)*.95));

raw.meetsThresh = raw.tProj > thresh;
clf
imagesc(raw.tProj, 'AlphaData', raw.meetsThresh)
colorbar
axis image;
title('thresholded tProjection')