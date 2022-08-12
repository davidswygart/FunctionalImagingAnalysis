path = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Ophthalmology\Research\SchwartzLab\Data\Zach Data\012722B';
lambda = [780, 825, 910, 980];
imgs = zeros(512,512,169,2*length(lambda));
for n=1:length(lambda)
[t,res{n},metadata{n},zs{n},ts{n},pos{n}] = fastLoadTiff(sprintf('%s%s%d_%05d.tif', path, filesep, lambda(n), 1));
%discard saturated pixels?
t = mean(t(:,:,1:2,:,:),5);
t = arrayfun(@(x) imgaussfilt3(t(:,:,x),res{n}), reshape(1:2*169,2,[]), 'UniformOutput', false);
imgs(:,:,:,2*n-1:2*n) = permute(reshape(cell2mat(t),[512,2,512,169]),[1,3,4,2]);
% imgs(:,:,:,2*n-1:2*n) = permute(imgaussfilt3(mean(t(:,:,1:2,:,:),5),res{n}),[1,2,4,3]);
end


imshow3D(cat(4,imgs(:,:,:,1:2),zeros(512,512,169,1))./1e4)


%SIFT keypoint detection

%DoG across locations and scales, to detect interest points invariance to
%scale and orientation

    % sequentially filter the image s+3 (?) times, increasing sigma by a factor
    % of 2^(1/s)
    
    % take the image that has twice the initial value of sigma (s+1?) and
    % downsample (standard anti-aliasing)
    
    % <- repeat
    
    % s = 3?

%determine location and scale at each candidate location

    % extrema in spatial and scale dimensions (27-n)

% select keypoints based on 'stability'

% assign oreintation(s) to each keypoint based on local gradient direction

% measure gradient at selected scale in region around each keypointw


for n=1:4
    pts{n} = detectSIFTFeatures(imgs(:,:,80,2*n-1)/1e4);
    [f{n},v{n}] = extractFeatures(imgs(:,:,80,2*n-1)/1e4,pts{n});

    subplot(3,4,n)
    imshow(imgs(:,:,80,2*n-1)/1e4)
    hold on
    plot(pts{n}.selectStrongest(25));

end
cnt=0;
for n=1:4
    for m=(n+1):4
        cnt = cnt + 1;
        ip{n,m} = matchFeatures(f{n},f{m});
        subplot(3,4,cnt+4)
        showMatchedFeatures(imgs(:,:,80,n),imgs(:,:,80,m),v{n}(ip{n,m}(:,1)),v{m}(ip{n,m}(:,2)))
    end
end

%align a range of slices in one volume to another
ref = 3;
plan = 153;
mov = 4;
set = 120:169;
scale = 1/1e4;
rows = 5;
cols = 4;

pts_ref = detectSIFTFeatures(imgs(:,:,plan,2*ref-1)*scale);
[f_ref,v_ref] = extractFeatures(imgs(:,:,plan,2*ref-1)*scale, pts_ref);
figure;
imshow(imgs(:,:,plan,2*ref-1)*scale)
hold on
plot(pts_ref.selectStrongest(100));

cnt = 0;
i = 0;
emp = cell(numel(set),1);
res = struct('features',emp,'pts',emp,'pairs',emp);
for s = set
    cnt = cnt+1;
    pts = detectSIFTFeatures(imgs(:,:,s,2*mov-1)*scale);
    [res(cnt).features,res(cnt).pts] = extractFeatures(imgs(:,:,s,2*mov-1)*scale,pts);
    res(cnt).pairs = matchFeatures(f_ref,res(cnt).features);
    
    if mod(cnt,rows*cols) == 1
        i = 1;
        figure;
    else
        i = i+1;
    end
    subplot(rows,cols,i);
    showMatchedFeatures(imgs(:,:,plan,2*ref-1)*scale, imgs(:,:,s,2*mov-1)*scale, v_ref(res(cnt).pairs(:,1)), res(cnt).pts(res(cnt).pairs(:,2)));
end






