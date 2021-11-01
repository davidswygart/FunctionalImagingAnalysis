function onsets = getTriggerTime(pixOrder, triggerImg)
linPixOrd = pixOrder(:);
linTrigImg(linPixOrd) = triggerImg(:);

t = linTrigImg > 1024; %no need to get fancy
onsets = find(t(2:end) & ~t(1:end-1));
%offsets = find(t(1:end-1) & ~t(2:end));

% figure
% plot(linImg)
% hold on
% scatter(onsets,linImg(onsets))

%move onset to end of line
% onsets = mod(onsets, pixPerLine) + 1 + onsets;
% offsets = mod(offsets, pixPerLine) + 1 + offsets;
end