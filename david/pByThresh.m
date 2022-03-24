function [p,v,tVals,pixPerc] = pByThresh(small,large,ThreshImg)
sortedThresh = sort(ThreshImg(:),'descend');

pixPerc = logspace(-3,0,50);

tInds = ceil(pixPerc*length(sortedThresh));
[tInds,isUnique] = unique(tInds, 'legacy'); %get rid of repeats if not a new pixel thresh
pixPerc = pixPerc(isUnique);
is3Pix = tInds>3; %Only run analysis if at least 3 pixels
tInds = tInds(is3Pix);
pixPerc = pixPerc(is3Pix);

tVals = sortedThresh(tInds);


p = nan(length(pixPerc),6);
v = nan(length(pixPerc),6);


for i = 1:length(pixPerc)
    % calculate si and shuffled si
    [siReal, siShuf] = analyzeSI(small, large, ThreshImg, tVals(i));

    % value of real
    v(i,1) = mad(siReal,0,1);
    v(i,2) = mad(siReal,1,1);
    v(i,3) = std(siReal,0,1);
    
    % p-value (mean MAD, med MAD, STD, Kurtosis)
    p(i,1) = calcP(mad(siReal,0,1), mad(siShuf,0,1), 'bigger');
    p(i,2) = calcP(mad(siReal,1,1), mad(siShuf,1,1), 'bigger');
    p(i,3) = calcP(std(siReal,0,1), std(siShuf,0,1), 'bigger');
    
    % calculate the difference between top and bottom 50%
    siReal = sort(siReal);
    siShuf = sort(siShuf,1);
    halfInd = floor(length(siReal)/2);
    
    smallHalfReal = siReal(1:halfInd);
    bigHalfReal = siReal(halfInd+1:end);
    
    smallHalfShuf = siShuf(1:halfInd,:);
    bigHalfShuf = siShuf(halfInd+1:end,:);
    
    avgDifReal = mean(bigHalfReal)-mean(smallHalfReal);
    avgDifShuf = mean(bigHalfShuf,1)-mean(smallHalfShuf,1);
    
    medDifReal = median(bigHalfReal)-median(smallHalfReal);
    medDifShuf = median(bigHalfShuf,1)-median(smallHalfShuf,1);
    
    v(i,4) = avgDifReal;
    v(i,5) = medDifReal;
    
    p(i,4) = calcP(avgDifReal, avgDifShuf, 'bigger');
    p(i,5) = calcP(medDifReal, medDifShuf, 'bigger');
    
    v(i,6) = mean(siReal);
    p(i,6) = calcP(mean(siReal), mean(siShuf,1), 'bigger');
    
    sprintf("completed %d of %d",i,length(pixPerc))
end

fsi = nan(6,1);
[r,c] = find(p<0.05);
for i = 1:6
    isSig = r(c==i);
    if isempty(isSig)
        fsi(i) = size(p,1);
    else
        fsi(i) = isSig(1);
    end
end

clf
subplot(2,1,2)
hold on
plot(pixPerc,p(:,1))
plot(pixPerc,p(:,2))
plot(pixPerc,p(:,3))
plot(pixPerc,p(:,4))
plot(pixPerc,p(:,5))
plot(pixPerc,p(:,6))

leg1 = sprintf('mean med (%.2f)', pixPerc(fsi(1)));
leg2 = sprintf('med med (%.2f)', pixPerc(fsi(2)));
leg3 = sprintf('std (%.2f)', pixPerc(fsi(3)));
leg4 = sprintf('half avg dif (%.2f)', pixPerc(fsi(4)));
leg5 = sprintf('half med dif (%.2f)', pixPerc(fsi(5)));
leg6 = sprintf('average si (%.2f)', pixPerc(fsi(6)));
legend(leg1,leg2,leg3,leg4,leg5,leg6,'Location','eastoutside')

ylabel('p value')
xlabel('pixel fraction')
set(gca, 'XScale', 'log')

subplot(2,1,1)
hold on
plot(pixPerc,v(:,1))
plot(pixPerc,v(:,2))
plot(pixPerc,v(:,3))
plot(pixPerc,v(:,4))
plot(pixPerc,v(:,5))
plot(pixPerc,v(:,6))

ylabel('measured value')
set(gca, 'XScale', 'log')

leg1 = sprintf('v = %.2f, t = %.2f', v(fsi(1),1),tVals(fsi(1)));
leg2 = sprintf('v = %.2f, t = %.2f', v(fsi(2),2),tVals(fsi(2)));
leg3 = sprintf('v = %.2f, t = %.2f', v(fsi(3),3),tVals(fsi(3)));
leg4 = sprintf('v = %.2f, t = %.2f', v(fsi(4),4),tVals(fsi(4)));
leg5 = sprintf('v = %.2f, t = %.2f', v(fsi(5),5),tVals(fsi(5)));
leg6 = sprintf('v = %.2f, t = %.2f', v(fsi(6),6),tVals(fsi(6)));
legend(leg1,leg2,leg3,leg4,leg5,leg6,'Location', 'eastoutside')

end