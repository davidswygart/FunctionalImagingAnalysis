%% I/O
region = 1;
input_path = 'E://MultiSMS/080421B/region%d_%05d.tif';
output_path = sprintf('E://MultiSMS/080421B/region%da', region);
cell_data = 'E:\MultiSMS\080421Bc4.mat';
h5_data = 'E:\MultiSMS\080421B.h5';
img = [];
ts=[];
saving = true;
spiking = true;
%region 2: [1,2] [3,4,5,6,7,8] [,9]
%TODO: group by sf, img size.... X.X
for i=2
    [timg, res, md, tts,pos] = fastLoadTiff(sprintf(input_path,region,i));
    [ny,nx] = size(timg,1,2);
    img = cat(3,img, permute(timg, [2, 1, 4, 3]));
    ts=cat(2,ts,tts);
end
sf = 1/mean(diff(ts));

%% Trigger
[on, off] = fastGetTriggerTime(img(:,:,:,end));
[~,e] = fastLoadSymphonyEpochAttributes(h5_data,'HexSMS');
if spiking
    s = matchSpikes(e,cell_data);
end
% [~,ei] = min(abs(ts(ceil(on./(nx*ny)))' - [e(:).t0]), [], 2);
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

% ti = nan(numel(onTime), max(preTime) + max(onTime) + 1);
% preTime(1) = 6;
% onTime(2) = 31;


ti = arrayfun(@(x,y,z)(-x:y) + z, ...
    preTime,...
    onTime,...
    ceil(on./(nx*ny))'+1,'uniformoutput',false);

spots = [[e(ei).offsetX] + [e(ei).cx]; [e(ei).offsetY] + [e(ei).cy]; e(ei).curSpotSize]';
[spots,~,ai] = unique(spots,'rows');

[sizes,~,w] = unique(spots(:,3));
N = accumarray(w,1); %number of spots per size

% epoched = arrayfun(@(x,y,z) cat(3,nan(ny-1,nx,mp-x,2),double(img(1:ny-1,:,z{1},[2,3])), nan(ny-1,nx,mo-y,2)),...
%     preTime,...
%     onTime,...
%     ti,...
%     'uniformoutput',false);
% epoched = permute(cat(5,epoched{:}),[1:3 5 4]);

%% Motion Correction
ops.useGPU = true;
ops.kriging = true;
ops.NiterPrealign = 20;
ops.Ly = ny-1; %ignore last row, due to projector artifact

xi = ceil(.175*nx) : floor(.825*nx);
ops.Lx = numel(xi);%nx; %to cut out artifact from pockel's cell


ops = alignIterative(double(img(1:ny-1,xi,:,2)),ops); %from Suite2P ~ align using transmitted IR?
reg = rigidRegFrames(double(img(1:ny-1,:,:,2)), ops, ops.dsprealign); %from Suite2P
regIR = rigidRegFrames(double(img(1:ny-1,:,:,3)), ops, ops.dsprealign); %from Suite2P

% ops = alignIterative(reshape(epoched(1:ny-1,xi,:,:,1),ny-1, ops.Lx, []),ops); %from Suite2P ~ align using transmitted IR?
% reg = rigidRegFrames(reshape(epoched(1:ny-1,:,:,:,1),ny-1, nx, []), ops, ops.dsprealign); %from Suite2P
% regIR = rigidRegFrames(reshape(epoched(1:ny-1,xi,:,:,2),ny-1, ops.Lx, []), ops, ops.dsprealign); %from Suite2P


%% ROI extraction
% ops = [];
% ops.useGPU = true;
% ops.fig = true;
% ops.nSVDforROI = inf;
% ops.diameter = 10; %!!!!
% % ops.yrange=res(2).*((1:ny-1) - (ny+1)/2);
% % ops.xrange=res(1).*((1:nx) - (nx+1)/2);
% ops.yrange = 1:ny-1; ops.xrange = 1:nx;
% ops.mimg1 = mean(img(1:ny-1,:,:,1),3);
% % ops.fs = sf;
% [stat,F,Fneu] = cellDetectionStandalone(img(1:ny-1,:,:,1), ops);
%
% epoched = reshape(F(:,ceil(on/(nx*ny)) + (-8:32)), size(F,1), length(on), []);
% dFoF = (epoched - mean(epoched(:,:,1:8),3)) ./ mean(epoched(:,:,1:8),3);
% epoched_neu = reshape(Fneu(:,round(on/(nx*ny)) + (-8:32)), size(F,1), length(on), []);
% dFoFneu = (epoched_neu - mean(epoched_neu(:,:,1:8),3)) ./ mean(epoched_neu(:,:,1:8),3);
% a = sum(dFoF(:,:,9:end).^2,3);

% epoched = double(reshape(reg(1:ny-1,:,fi + window), nx*(ny-1), length(fi), length(window)));
% epoched = double(reshape(img(1:ny-1,:,ceil(on/(nx*ny)) + window,1), nx*(ny-1), length(on), length(window)));
% dFoF = (epoched - mean(epoched(:,:,window<0),3)) ./ mean(epoched(:,:,window<0),3);
% epoched = reshape(double(img(1:ny-1,:,ti,1)), [(ny-1), nx, size(ti)]);

% epoched = arrayfun(@(x,y,z) cat(3,nan(ny-1,nx,mp-x),double(reg(1:ny-1,:,z{1})), nan(ny-1,nx,mo-y)),...
%     preTime,...
%     onTime,...
%     ti,...
%     'uniformoutput',false);
% epoched = cat(4,epoched{:});
%
% bl = nanmean(epoched(:,:,1:mp,:),3);
% dFoF = (epoched - bl) ./ (abs(bl) + eps); %this is because the noise is 0 mean...
% a = sum(dFoF(:,:,mp+1:end,:).^2,3);
% at = splitapply(@(x) mean(x,2), reshape(a,(ny-1)*nx,[]), ai'); %validate this...

% at(57,:)= []; %TODO: no examples of this one, so it shows up NaN?
%

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
% sms = zeros(12,1);
% % sta = zeros(12,2);
% for j = 1:size(dFoF,1)
%     at = accumarray(ai,a(j,:)',[],@mean);
%     at = at - min(at); at = at./max(at);
%
% end

if saving
    try
%         rmdir('E:/MultiSMS/080421B/region3a','s')
        rmdir(output_path,'s');
    end
%     mkdir('E:/MultiSMS/080421B/region3a')
%     mkdir('E:/MultiSMS/080421B/region3a/rois')
%     mkdir('E:/MultiSMS/080421B/region3a/rois/traces')
    mkdir(output_path)
    mkdir(sprintf('%s/rois',output_path))
    mkdir(sprintf('%s/rois/traces',output_path))
    
end
%%
if saving
    sta = zeros(12,2,2);
    figure('units','normalized','position',[0,0,1,1]);clf;set(gcf,'color','w');
    for roi = 1:nROIs
        clf;
        loc = res.*(stat(roi).ellipse-[ny,nx]/2-.5);
        for o=1:2
            if o == 1
                at = stimOnAvg;
            else
                at = stimOffAvg;
            end
            for i = 1:numel(sizes)
                si = spots(:,3)==sizes(i);
                subplot(6,4,i+(o-1)*12);
                caxis([0,max(at(:,roi))]);
                hexPatch(spots(si,:), abs(at(ai(si),roi)));
                caxis([0,max(at(:,roi))]);
                axis equal;
                
                att = at(ai(si),roi);
                att(isnan(att)) = 0;
                
                if sum(att)>0
                    sta(i,o,:) = att'*spots(si,1:2)./sum(att);% * sqrt(N(i));
                elseif nnz(si)>1
                    %         sta(i,:) = mean(spots(si,1:2),1);
                else
                    sta(i,o,:) = spots(si,1:2);
                end
                hold on;
                plot(sta(i,o,1),sta(i,o,2),'r+');
                %     loc = res.*(stat(roi).med-[ny,nx]/2-.5);
                plot(loc(:,2),-loc(:,1),'g');
                %     plot(mu(1),mu(2),'g+');
            end
        end
        print(sprintf('%s/rois/roi_%02d',output_path,roi),'-dpng');
    end
end
%%
if saving
    figure('units','normalized','position',[0,0,1,1]);clf;set(gcf,'color','w');
    spots = [[e(ei).offsetX] + [e(ei).cx]; [e(ei).offsetY] + [e(ei).cy]; e(ei).curSpotSize]';
    [ai, x, y, spotSize] = findgroups(spots(:,1), spots(:,2), spots(:,3));
    % [spots,~,ai] = unique(spots,'rows');
    for roi = 1:nROIs
        clf;
        loc = res.*(stat(roi).ellipse-[ny,nx]/2-.5);
        %     mkdir(sprintf('E:/MultiSMS/080421B/rois/roi_%02d',roi));
        att = squeeze(dFoF(roi,:,:) ./ max(abs(dFoF(roi,:,:)),[],'all'))';
        for i = 1:numel(sizes)
            si = spotSize==sizes(i);
            ii = ismember(ai,find(si));
            
            %     subplot(3,4,i);
            cla;
            %     hold on;
            
            hexTrace(spots(ii,:), att(ii,:), [mp mp+sf]); %???
            axis off;
            hold on;
            plot(loc(:,2),-loc(:,1),'g');
            title(sprintf('ROI %02d, Spot diameter = %.2f\\mum',roi,sizes(i)))
            print(sprintf('%s/rois/traces/roi_%02d_size_%02d',output_path,roi,i),'-dpng');
            
        end
    end
end
%%
if spiking && saving
    
    mkdir(sprintf('%s/rasters',output_path));
    roi = 7; %index of the ROI to associate with spikes
    figure('units','normalized','position',[0,0,1,1]);clf;set(gcf,'color','w');
    spots = [[e(ei).offsetX] + [e(ei).cx]; [e(ei).offsetY] + [e(ei).cy]; e(ei).curSpotSize]';
    [ai, x, y, spotSize] = findgroups(spots(:,1), spots(:,2), spots(:,3));
    loc = res.*(stat(roi).ellipse-[ny,nx]/2-.5);
    att = squeeze(dFoF(roi,:,:) ./ max(abs(dFoF(roi,:,:)),[],'all'))';
    for i = 1:numel(sizes)
        si = spotSize==sizes(i);
        ii = ismember(ai,find(si));
        
        %     subplot(3,4,i);
        cla;
        %     hold on;
        
        hexTrace(spots(ii,:), att(ii,:), [mp mp+sf], spikes(ii), [500 1000]); %???
        axis off;
        hold on;
        plot(loc(:,2),-loc(:,1),'g');
        title(sprintf('ROI %02d, Spot diameter = %.2f\\mum',roi,sizes(i)))
        print(sprintf('%s/rasters/roi_%02d_size_%02d',output_path,roi,i),'-dpng');
        
    end
end
%%
% end
%
% mu_hat = log(N)'*sta / sum(log(N)); %weight the grids by information content
% for i = 1:12
%     si = spots(:,3)==sizes(i);
%     if nnz(si)>1
%         f = scatteredInterpolant(spots(si,1), spots(si,2), a(si),'natural');
%         sms(i) = f(mu(1),mu(2));
%     else
%         sms(i) = a(si);
%     end
% end
% figure(2);clf;
% % semilogx(sizes, trueSMS); hold on; semilogx(sizes,sms./sms(1)); scatter(spots(:,3),a./sms(1),'+','markeredgecolor',[.7,.7,.7]);
% figure(1);

%% Notes
%ROI detection

%goals:
%remove neuropil signal -- spatially smooth, basis functions
%signal at pixel k: sum(B(k,j)*n(j,t), j)
%isotropic, 2d raised cosine, default spacing ~ 3x cell diameter
%isolate signal from each roi i, f(i,t)
%each roi pixel k scaled by non-negative constant, lambda(k,i)

%full model for recorded signal r(k,t):
%   r(k,t) = sum(lambda(k,i)*f(i,t), i) + sum(B(k,j)*n(j,t),j) + eta(k,t)
%   where eta ~ N(0,s^2)

%approach:
%   1) reduce to ~2k components (from 4k) using svd
%   2) detect new ROIs
%   3) estimate ROI/neuropil timecourses and residuals
%   4) estimate spatial distribution of ROIs
%   5) return to 2, iterate 10-30x?

%% initialize model
%lambda = sparse 0
%B = cosine functions...
%f = zeros, roi signals
%n = zeros, neuropil signals
%kernel_w = ?
%mindist = 2*avg_soma_size_in_px
%roi_count = 0

%% main loop
%for iter:
%   %% Source detection
%   resid = raw - model(lambda,B,f,n)
%   map = fastCorrelationMap(resid,kernel_w)
%   new_rois = findpeaks_descending(map, mindist)
%   for i in new_rois:
%       lambda(new_rois(i),i+roi_count) = map(new_rois(i))
%   roi_count+=new_roi_count
%
%   %% Source extraction
%   A = [lambda, zeros; zeros, B]
%   fn = pseudoinv(A)*raw
%   [f,n] = [fn(1:F), fb(F+1:end)]
%
%   %% Source re-assignment
%   signals = raw - B*n
%   for i in roi_count:
%       others = setdiff(1:roi_count, i)
%       signal = signals - lambda(:,others)*f(others,:)
%       for pixels k:
%           lambda(k,i) = lsq(x=f(i,:), y= signal(k,:) )
