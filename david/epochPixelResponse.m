function r = epochPixelResponse(green, time)
%calculate the response for each pixel for each epoch

%% indicate pre and stim times
isPre = time < 0;
isStim = time > 0 & time < 1;

%% calculate pre-time values for each epoch
temp = green;
temp(~isPre) = nan;
r.pre = squeeze(mean(temp, 3, 'omitnan'));
r.stdPre = squeeze(std(temp,0, 3, 'omitnan'));

%% calculate stim-time values for each epoch
%for each pixel
temp = green;
temp(~isStim) = nan;
r.stim = squeeze(mean(temp, 3, 'omitnan'));
r.stdStim = squeeze(std(temp,0, 3, 'omitnan'));

%% calculate pixel-wise measures of response
r.dF = r.stim - r.pre;
r.dPrime = r.dF ./ ((r.stdStim + r.stdPre) / 2);
r.dFoF = r.dF ./ r.pre;
end