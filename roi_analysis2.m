%% I/O
input_path = 'E://MultiSMS/080421B/region1_00002.tif';
h5_data = 'E:\MultiSMS\080421B.h5';
output_path = 'E:\MultiSMS\080421B\region1\manual_rois_2\';
rois_path = 'E:\MultiSMS\080421B\region1_00002_rois.zip'; 

img = [];
stim = [];
ts = [];
for fo = dir(input_path)'
    [timg, res, ~, ~, tts, ~] = fastLoadTiff(fullfile(fo.folder, fo.name));
    [ny,nx] = size(timg,1,2);
    stim = cat(3,stim,squeeze(timg(:,:,end,:)));
    ts = cat(1, ts, tts);
% end
% img = permute(img,[2,1,4,3]);
% 
sf = 1/median(diff(ts));
    %% Resample the channels and save as .h5
%     out_file = sprintf('%s%s.h5',output_path, fo.name);
%     h5create(out_file, '/data', size(timg,[1, 1, 3, 4]), 'Datatype','uint16','ChunkSize', [size(timg,[1,1,3]),200]);
% 
%     for n=1:size(timg,3)
%      
%         F = griddedInterpolant({1:ny, 1:nx}, double(squeeze(timg(:,:,n,:))));
%         F.Method = 'nearest';
%         F.ExtrapolationMethod = 'nearest';
% 
%         rs = F({(1:ny)', linspace(.5,nx+.5,ny)'});
%         h5write(out_file, '/data', reshape(rs, ny, ny, 1, []), [1, 1, n, 1], [ny, ny, 1,  size(timg,4)]);
%     end
    F = griddedInterpolant({1:ny, 1:nx}, double(squeeze(timg(:,:,1,:))));
    F.Method = 'nearest';
    F.ExtrapolationMethod = 'nearest';

    img = cat(3,img,F({(1:ny)', linspace(.5,nx+.5,ny)'}));
end
%%
ops = [];
ops.useGPU = true;
ops.kriging = true;
ops.NiterPrealign = 20;
ops.Ly = ny-1; %ignore last row, due to projector artifact
% 
xi = ceil(.175*ny) : floor(.825*ny);
ops.Lx = numel(xi);%nx; %to cut out artifact from pockel's cell
% 
ops = alignIterative(double(img(1:ny-1,xi,:,1)),ops); %from Suite2P ~ align using transmitted IR?
reg = rigidRegFrames(double(img(1:ny-1,:,:,1)), ops, ops.dsprealign); %from Suite2P
% regIR = rigidRegFrames(double(img(1:ny-1,:,1:4e4,3)), ops, ops.dsprealign); %from Suite2P
save(sprintf('%sreg1.mat',output_path),'reg','-v7.3');

%% Post- motion correction...
% img = fastLoadTiffGen('E:\MultiSMS\042722B\resampled\chirp\suite2p\plane0\reg_tif\file000_chan0.tif');
% [on, off] = fastGetTriggerTime(stim); %subframe timing
[~,e] = fastLoadSymphonyEpochAttributes(h5_data,'HexSMS');

trigTime = posixtime(datetime([e(:).t0],'convertfrom','posixtime') + milliseconds([e(:).preTime]))';
[ont,on] = min(abs(ts-trigTime'),[],1);
oni = ont<1;
e = e(oni);
on = on(oni);


mp = floor(e(1).preTime*1e-3 * sf);
mo = ceil((e(1).stimTime)*1e-3 * sf);
ei = -mp:ceil((e(1).stimTime + e(1).tailTime)*1e-3 * sf); %assumes the same epoching...
ti = on + ei';

ROIs= readIJROIs(rois_path);

trend_time = ceil(sf*20);
% fn = @(x) padarray(mean(x,1),[0,trend_time],
Fe = zeros(numel(ROIs), numel(ei), numel(on));
F = zeros(size(img, 3), numel(ROIs));
for n=1:numel(ROIs)
    %applies a function to the data (pix -by- time)
    F(:,n) = ROIs(n).apply(@mean,reg);
%     pF = padarray(F(:,n), [trend_time,0], F(1,n),'pre');
%     trend = arrayfun(@(x) prctile(pF(x-trend_time:x),8), trend_time+1:length(pF));
%     Fi = (F(:,n) - trend') ./ trend';
%     Fe(n,:,:) = reshape(Fi(ti), size(ti));
    Fe(n,:,:) = reshape(F(ti,n), size(ti));
end
bl = mean(Fe(:,1:mp,:),2);
dFoF = (Fe - bl) ./ bl;

spots = [[e.offsetX] + [e.cx]; [e.offsetY] + [e.cy]; e.curSpotSize]';
[ai, x, y, spotSize] = findgroups(spots(:,1), spots(:,2), spots(:,3));

[sizes,~,w] = unique(spots(:,3));


n = 16;
i = 1;
for n=1:numel(ROIs)
att = squeeze(dFoF(n,:,:) ./ max(abs(dFoF(n,:,:)),[],'all'))';

for i=1:12
si = spotSize==sizes(i);
ii = ismember(ai,find(si));

%     subplot(3,4,i);
cla;
%     hold on;

hexTrace2(spots(ii,:), att(ii,:), [mp mp+sf]); %???
axis off;
hold on;
%             plot(loc(:,2),-loc(:,1),'g');
title(sprintf('ROI %02d, Spot diameter = %.2f\\mum',n,sizes(i)),'color','w')
    print(sprintf('%s/roi_%02d_size_%02d',output_path,n,i),'-dpng');
end
end