function [onsets,offsets] = fastGetTriggerTime(image)
%this assumes that the image has been permuted such that it is in column
%major order....
%please keep in mind that sub-line precision is probably false precision

%first convert from int16 to uint16 to truncate negative values to 0
%then right shift by 13 to remove noise
% timing = reshape(permute(uint8(bitshift(uint16(image),-13)),[2,1,3]),[],1);
%we convert to a temporal signal by permuting back to temporal order and
%linearizing

% onoff = [find(timing>0 & timing([1, 1:end-1])==0),find(timing>0 & timing([2:end, end])==0)];

%above is a bit slower than the naive method:
timing = reshape(permute(image,[2,1,3]),[],1);

% onoff = [...
%     find(timing>1024 & timing([1, 1:end-1])<=1024),...
%     find(timing>1024 & timing([2:end, end])<=1024)...
%     ];

t = timing > 1024; %no need to get fancy
onoff = [find(t(2:end) & ~t(1:end-1)), find(t(1:end-1) & ~t(2:end))];
%doing the ~ operator twice is faster than saving in a separate variable,
%but not by much

keep = onoff(:,2) - onoff(:,1) > 128; %debouncing
onsets = onoff(keep,1);
offsets = onoff(keep,2);



end