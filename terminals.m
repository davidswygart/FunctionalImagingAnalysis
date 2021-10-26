%% I/O
input_path = 'D://Images/102121B/102121Bc2/102121Bc2_SMS.tif';
cell_data = 'D://Images/102121B/102121Bc2.mat';
h5_data = 'D://Images/102121B/102121B.h5';
img = [];
ts=[];
spiking = true;
for i=1
    [timg, res, md, tts,pos] = fastLoadTiff(sprintf(input_path,i));
    [ny,nx] = size(timg,1,2);
    img = cat(3,img, permute(timg, [2, 1, 4, 3]));
    ts=cat(2,ts,tts);
end
sf = 1/mean(diff(ts));

%% Trigger
[on, off] = fastGetTriggerTime(img(:,:,:,end));
[~,e] = fastLoadSymphonyEpochAttributes(h5_data,'SpotsMultiSize');
if spiking
    s = matchSpikes(e,cell_data);
end
trigTime = posixtime(datetime([e(:).t0],'convertfrom','posixtime') + milliseconds([e(:).preTime]))';
delta = ts(ceil(on./(nx*ny))+1) - trigTime;
delta(abs(delta) > 1/sf) = inf;
rm = all(isinf(delta),1);
[~,ei] = min(delta,[],1);
ei(rm) = [];
on(rm) = [];

preTime = floor([e(ei).preTime]*1e-3 * sf);
onTime = floor(([e(ei).stimTime]+[e(ei).tailTime])*1e-3 * sf);
mp = max(preTime);
mo = max(onTime);

if spiking
    spikes = s(ei);
end


ti = arrayfun(@(x,y,z)(-x:y) + z, ...
    preTime,...
    onTime,...
    ceil(on./(nx*ny))'+1,'uniformoutput',false);

spots = [e(ei).curSpotSize]';
[spots,~,ai] = unique(spots,'rows');

[sizes,~,w] = unique(spots(:,end));
N = accumarray(w,1); %number of spots per size

%% Motion Correction
ops.useGPU = true;
ops.kriging = true;
ops.NiterPrealign = 20;
ops.Ly = ny-1; %ignore last row, due to projector artifact

xi = ceil(.175*nx) : floor(.825*nx);
ops.Lx = numel(xi);%nx; %to cut out artifact from pockel's cell


ops = alignIterative(double(img(1:ny-1,xi,:,3)),ops); %from Suite2P ~ align using transmitted IR?
reg = rigidRegFrames(double(img(1:ny-1,:,:,1)), ops, ops.dsprealign); %from Suite2P
regIR = rigidRegFrames(double(img(1:ny-1,:,:,3)), ops, ops.dsprealign); %from Suite2P

% ops = alignIterative(reshape(epoched(1:ny-1,xi,:,:,1),ny-1, ops.Lx, []),ops); %from Suite2P ~ align using transmitted IR?
% reg = rigidRegFrames(reshape(epoched(1:ny-1,:,:,:,1),ny-1, nx, []), ops, ops.dsprealign); %from Suite2P
% regIR = rigidRegFrames(reshape(epoched(1:ny-1,xi,:,:,2),ny-1, ops.Lx, []), ops, ops.dsprealign); %from Suite2P


%% dFoF model?
ops = [];
ops.useGPU = true;
ops.fig = true;
ops.nSVDforROI = inf;
ops.diameter = 10; %!!!!
ops.yrange = 1:ny-1; ops.xrange = 1:nx;
ops.fs = sf;
% ops.mimg1 = mean(reshape(epoched(1:ny-1,:,mp+1:min(onTime),:),ny-1, nx, []),3);

ops.mimg1 = mean(reg,3);
% [stat,F,Fneu] = cellDetectionStandalone(reshape(epoched(1:ny-1,:,mp+1:min(onTime),:),ny-1, nx, []), ops);

[stat,F,Fneu] = cellDetectionStandalone(reg, ops);
nROIs = size(F,1);
epoched = arrayfun(@(x,y,z) cat(2,nan(nROIs,mp-x),double(F(:,z{1})), nan(nROIs,mo-y)),...
    preTime,...
    onTime,...
    ti,...
    'uniformoutput',false);
epoched = cat(3,epoched{:}) - min(F,[],'all');

bl = nanmean(epoched(:,1:mp,:),2);
dFoF = (epoched - bl) ./ (abs(bl) + eps); %this is because the noise is 0 mean...

a = splitapply(@(x) mean(x,2), shiftdim(dFoF,1), ai');

stimOn = sum(dFoF(:,mp+(1:floor(sf)),:).^2,2);
stimOnAvg = splitapply(@(x) mean(x,2), reshape(stimOn,nROIs,[]), ai')'; %validate this...

stimOff = sum(dFoF(:,(mp+ceil(sf)):end,:).^2,2);
stimOffAvg = splitapply(@(x) mean(x,2), reshape(stimOff,nROIs,[]), ai')'; %validate this...


%% Response mapping
