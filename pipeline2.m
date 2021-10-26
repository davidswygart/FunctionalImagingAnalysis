%% I/O
file_name_and_path = 'E://MultiSMS/061821B/region2_cell3_00004.tif';
[img, res] = fastLoadTiff(file_name_and_path);
[ny,nx] = size(img,1,2);
img = permute(img, [2, 1, 4, 3]);

%% Bidirectional scan phase correction
%Suite2P has a function for this, but seems like scanimage did well enough

%% Image registration -- slow
ops.useGPU = true;
ops.kriging = true;
ops.NiterPrealign = 20;
ops.Lx = ny;
ops.Ly = nx-1; %ignore last row, due to projector artifact
ops = alignIterative(double(img(1:nx-1,22:nx-23,:,3)),ops); %from Suite2P ~ align using transmitted IR?
reg = rigidRegFrames(double(img(1:nx-1,:,:,2)), ops, ops.dsprealign); %from Suite2P
%NOTE: not clear if it's better to align from the green or IR channel
% IR has higher signal, but worse optical sectioning

%% Deconvolution -- slow ~~ before/after registration? or not at all?
psf = fspecial('gaussian',12,0.9); %just a guess
reg_dc = deconvlucy( reg, psf);

%% ROI detection
map = fastCorrelationMap(reg_dc,0.9); %adapted from Suite2P paper
%blur parameter was chosen arbirarily

%% Epoch windowing
[on, off] = fastGetTriggerTime(img(:,:,:,end));
off(end) = []; on(end) = []; %there's a light step at the end to kill the mean?
% dur = floor(min(off-on)/(nx*ny));
all_E = reshape(reg_dc(:,:,round(on/(nx*ny)) + (-8:32)), nx, ny, length(on), []);
all_dFoF = (all_E - mean(all_E(:,:,:,1:8), 4)) ./ mean(all_E(:,:,:,1:8),4);


%% Suppression index stuff
window = [-30 30]; %pretty arbitrary
[on,off] = fastGetTriggerTime(img(:,:,:,3));
all_E = reshape(reg_dc(:,:,round(on/(nx*ny)) + (window(1):window(2))), nx-1, ny, length(on),[]);
all_dFoF = (all_E - mean(all_E(:,:,:,1:-window(1)),4)) ./ mean(all_E(:,:,:,1:-window(1)),4);

small = mean(all_dFoF(:,:,1:2:end - mod(length(on),2),1-window(1):end).^2,4);
large = mean(all_dFoF(:,:,2:2:end,1-window(1):end).^2,4);

% NOTE: intentionally comparing subsequent epochs to each other, to curb
% effect of slow drifts in baseline fluorescence, bleaching
all_SI = (small - large) ./ (small + large);

all_SI_map = squeeze(median(all_SI,3));


%% Plot SI vs ROI
roi_thresh = .1; %crude, but a simple start
% note that as written this will bias us away from detecting terminals with
% SI variance

masked = nan(nx-1,ny);
masked(map > roi_thresh) = all_SI_map(map > roi_thresh);
% imagesc(masked)
