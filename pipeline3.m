%% I/O
[img, res,md,zs,ts,pos] = fastLoadTiff('E://MultiSMS/081821B/region0_cell1_00002.tif');
[ny,nx] = size(img,1,2);
[on, off] = fastGetTriggerTime(squeeze(img(:,:,end,:)));
sf = 1./median(diff(ts));


ROIs = readIJROIs('E:\MultiSMS\081821B\region0_cell1_00002.roi');

[~,epochs] = fastLoadSymphonyEpochAttributes('E://MultiSMS/081821B.h5','HexSMS');
% TODO: load from db....

g_raw = squeeze(img(1:ny-1,1:nx,2,:));
%% segmentation
seg = zeros(numel(ROIs),size(g_raw,3));
[qx,qy] = meshgrid(1:ny-1,1:nx);
for n=1:numel(ROIs)
    seg(n,:) = ROIs(n).apply(@mean, g_raw,qx,qy);
end


%% epoching
%TODO: pretime, stimtime, tail time...
window = floor(-sf/2) : ceil(2*sf);
t0 = [nnz(window<=0), nnz(window<=0) + ceil(sf)];
onset_frame = floor(on/(nx*ny));
%TODO: with line clocking, we want to round up or down differently for
%lines above and below...

epoched = reshape(seg(:,onset_frame + window),[numel(ROIs),numel(on), numel(window)]);

[~, ei] = min(abs([epochs(:).t0] - ts(onset_frame)'),[],2);

params = [epochs(ei).cx ; epochs(ei).cy; epochs(ei).curSpotSize]';
[~,~,g] = unique(params,'rows');
[sizes,~,si] = unique(params(:,3));


bl = mean(epoched(:,:,window<0),3);
dFoF = (epoched - bl) ./ bl;


close all;
for m=1:numel(sizes)
    figure('windowstyle','docked');
    set(gcf,'color','w');
end
for n=1:numel(ROIs)
for m=1:numel(sizes)
figure(m);clf;axis off;
hold on
fill(-nx/2 + ROIs(n).coords(:,1),nx/2-ROIs(n).coords(:,2),'g')
hexTrace(params(si==m,:), squeeze(dFoF(n,si==m,:)./max(abs(dFoF(n,si==m,:)),[],'all')),t0)
% saveas(gcf, sprintf('E:\\MultiSMS\\080421B\\region1\\manual_rois\\roi_%02d_size_%02d.png',n,m));
end
end

