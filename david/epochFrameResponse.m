function r = epochFrameResponse(green, time, thresh)
%calculate the response everything above threshold

%% make threshold logical size of image

[y,x,t,e] = size(green);
thresh = repmat(thresh,[1,1,t,e]);


%% indicate pre and stim times
isPre = time < 0;
isStim = time > 0 & time < 1;

%% calculate pre-time values for each epoch
temp = green;
temp(~isPre | ~thresh) = nan;
r.pre = squeeze(mean(temp, [1,2,3], 'omitnan')); %average accross space and time, leaving only epochs
r.stdPre = squeeze(std(temp,0, [1,2,3], 'omitnan'));

%% calculate stim-time values for each epoch (and uncertainty)
temp = green;
temp(~isStim | ~thresh) = nan;
r.stim = squeeze(mean(temp, [1,2,3], 'omitnan'));
r.stdStim = squeeze(std(temp,0, [1,2,3], 'omitnan'));

%% calculate dF
r.dF = r.stim - r.pre;
c = cov(r.stim, r.pre);
r.df_std = sqrt(r.stdStim.^2 + r.stdPre.^2 - 2*c(2,1)); %propogation of error

%% calculate dPrime
r.dPrime = r.dF ./ ((r.stdStim + r.stdPre) / 2);

%% calculate dFoF
r.dFoF = r.dF ./ r.pre;
c = cov(r.dF, r.pre);

r.dFoF_std = abs(r.dFoF) .* sqrt( (r.df_std./r.dF).^2  + (r.pre./r.stdPre).^2 - (2 .* c(1,2) ./ r.dF ./ r.pre)  ) ;%propogation of error
r.dFoF_sem = r.dFoF_std ./ sqrt(squeeze(sum(thresh, [1,2,3])));
end