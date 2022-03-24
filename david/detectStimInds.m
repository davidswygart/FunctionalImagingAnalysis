function [onInd, offInd] = detectStimInds(stim, isBi)
    %if obj.SI.hScan2D.bidirectional
    %end
    stim = double(stim);
    [nR,nC,nZ] = size(stim);
    
    npixInFrame = nR*nC;
    
    if isBi
        stim(2:2:end,:,:) = flip(stim(2:2:end,:,:), 2);
    end
    
    stim = permute(stim, [2,1,3]);
    linImg = reshape(stim, [],1);
    
    linImg = linImg - min(linImg);
    thresh = max(linImg) / 2;
    
    linImg(linImg<(thresh)) = 0;
    linImg(linImg>=(thresh)) = 1;
    
    
    dif = diff(linImg);

    on = find(dif == 1);
    off = find(dif == -1);
    
    %get rid of any on/off detections within 8 frames of the image end
    framesFromEnd = (length(linImg)-off) / npixInFrame;
    off = off(framesFromEnd > 8);
    framesFromEnd = (length(linImg)-on) / npixInFrame;
    on = on(framesFromEnd > 8);

    
    on = on(linImg(on+npixInFrame) == 1);  % only include ons for which the next frame is actually on
    off = off(linImg(off+npixInFrame) == 0);% only include offs for which the next frame is actually off
    
    on(diff(on) < npixInFrame) = [];% throw away ones that occur within 1 frame of the previous detection
    off(diff(off) < npixInFrame) = [];% throw away ones that occur within 1 frame of the previous detection
    

    
    %% Error checking
    if length(on) > length(off)
        warning('More Onsets were detected than Offsets, removing last onset')
        on = on(1:end-1);
    end

    if min(diff(on)) < npixInFrame*10
        warning('The the minimum time between detected epochs is less than a 10 frames')
    end
    
    if min(off - on) < 0
        warning('An epoch offset is detected before the epoch onset')
    end
    
    if (min(off - on) > 0) && (min(off - on) < npixInFrame)
        warning('The difference between onset and offset of an epoch is less than a single frame')
    end
    
%     hold on
%     plot(linImg,'k')
%     scatter(on,linImg(on))
%     scatter(off,linImg(off))
    
    %%
    [c,r,z] = ind2sub([nC, nR, nZ], on);
    %shift down to next line because projector only turns on during
    %blanking period
    for i = 1:length(r)
        if r(i) < nR
            r(i) = r(i) + 1;
        else
            r(i) = 1;
            z(i) = z(i)+1;
        end
    end
        
    %choose start of line ()
    oddRows = logical(bitget(abs(r),1));
    c(oddRows) = 1;
    c(~oddRows) = nC;
            
    onInd = [r,c,z];
    
    
    [c,r,z] = ind2sub([nC, nR, nZ], off);
        %shift down to next line because projector only turns on during
    %blanking period
    for i = 1:length(r)
        if r(i) < nR
            r(i) = r(i) + 1;
        else
            r(i) = 1;
            z(i) = z(i)+1;
        end
    end
        
    %choose start of line ()
    oddRows = logical(bitget(abs(r),1));
    c(oddRows) = 1;
    c(~oddRows) = nC;
    offInd = [r,c,z];
    
end