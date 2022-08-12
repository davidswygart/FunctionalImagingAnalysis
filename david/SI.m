%% calculate SI using best spot
big = spot(end).avgDFoF;
small = spot(bestInd).avgDFoF;

% big_err = spot(end).stdDFoF;
% small_err = spot(bestInd).stdDFoF;
% 
% num = (small - big);
% dem = (small + big);
% 
% num_err = small_err.^2 + big_err.^2;
% den_err = num_err;

si = (small - big) ./ (small + big);


supp = struct;
supp.avgSI = mean(si,3);
supp.stdSI = std(si,0,3);
% supp.stdSI = errMult(num,dem,num_err,den_err);

clf
%imagesc(supp.avgSI)
%imagesc(supp.avgSI, 'AlphaData', raw.meetsThresh)
imagesc(supp.stdsi)
caxis([-1 2])
colorbar

%%

big_err = spot.pix.dFoF_err(:,:,end);

for i = 1:spot.n
%     A = spot.pix.avgStim(:,:,i);
%     a = spot.pix.avgPre(:,:,i);  
%     sA = spot.pix.stdStim(:,:,i);
%     sa = spot.pix.stdPre(:,:,i); 
%     
%     [si, err] = calcSI_DFoF(A,a,B,b,sA,sa,sB,sb);
%     %[si, err] = calcSI_DF(A,a,B,b,sA,sa,sB,sb);
%     
%     spot.pix.dFoF
%     
%     spot.pix.avgSI(:,:,i) = si;
%     spot.pix.stdSI(:,:,i) = err;
    
    small = spot.pix.dFoF(:,:,i);
    small_err = spot.pix.dFoF_err(:,:,i);
    
    num = (small - big);
    dem = (small + big);
    spot.pix.avgSI(:,:,i) = num ./ dem;
    
    num_err = small_err.^2 + big_err.^2;
    den_err = num_err;
    
    spot.pix.stdSI(:,:,i) = errMult(num,dem,num_err,den_err);
end



%% display SI for every size
clf
im3d = spot.pix.avgSI;
%im3d = spot.pix.stdSI;

%im3d = im3d .* repmat(raw.meetsThresh ,1,1,spot.n); % apply a threshold

imshow3D(im3d, [-1.5 2])
colorbar
axis image;

%% show SI for the best spot
clf

spotInd = bestSpotInd
%im2D = spot.pix.stdSI(:,:,spotInd);
im2D = spot.pix.avgSI(:,:,spotInd);

%imagesc(im2D, 'AlphaData', raw.meetsThresh)
imagesc(im2D)
caxis([-1 1])
axis image;
colorbar

%% histogram of the above values only at thresholded pixels
clf
edges = linspace(0,2,20);
histogram(im2D(raw.meetsThresh), edges)

%% histogram of dPrime values between pixel SI
err = spot.pix.stdSI(:,:,spotInd);
err = err(raw.meetsThresh);
avg = spot.pix.avgSI(:,:,spotInd);
avg = avg(raw.meetsThresh);

df = abs(avg - avg');

avgErr = (err + err') / 2;

dPrime = df ./ avgErr;

isUnique = logical(tril(ones(size(dPrime)),-1));
uniqueDprime = dPrime(isUnique);

edges = linspace(0,30,20);
histogram(uniqueDprime, edges)
