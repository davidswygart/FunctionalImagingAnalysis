%% I/O
file_name_and_path = '../031721B/031721_00010.tif';
[img, res] = fastLoadTiff(file_name_and_path);
[ny,nx] = size(img,1,2);
img = permute(reshape(img,ny, nx, 3,[]), [2, 1, 4, 3]);

%% Bidirectional scan phase correction
%Suite2P has a function for this, but seems like scanimage did well enough

%% Image registration
ops.useGPU = true;
ops.kriging = true;
ops.NiterPrealign = 20;
ops.Lx = ny;
ops.Ly = nx-1; %ignore last row, due to projector artifact
ops = alignIterative(double(img(1:nx-1,:,:,2)),ops); %from Suite2P ~ align using transmitted IR?
reg = rigidRegFrames(double(img(1:nx-1,:,:,1)), ops, ops.dsprealign); %from Suite2P
%NOTE: not clear if it's better to align from the green or IR channel
% IR has higher signal, but worse optical sectioning

%% Deconvolution ~~ before/after registration? or not at all?
psf = fspecial('gaussian',12,0.9); %just a guess
reg_dc = deconvlucy( reg, psf);

%% ROI detection
map = fastCorrelationMap(reg_dc,0.9); %adapted from Suite2P paper

%% Epoch windowing
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
masked = nan(nx-1,ny);
masked(map > roi_thresh) = all_SI_map(map > roi_thresh);
% imagesc(masked)
