function [p,difference,tVals] = diffByThresh(small,large,ThreshImg)
ThreshImg(isnan(ThreshImg)) = -inf;
sortedThresh = sort(ThreshImg(:),'descend');

tVals = 1 - [0:0.1:1];
pixCounts = sum(sortedThresh > tVals,1);
tVals(pixCounts < 4) = []; % get rid of threshold values with less than 4 pixels



topAvg = nan(length(tVals),1);
bottomAvg = nan(length(tVals),1);
difference = nan(length(tVals),1);
p = nan(length(tVals),1);


for i = 1:length(tVals)
    % calculate si and shuffled si
    [siReal, siShuf] = analyzeSI(small, large, ThreshImg, tVals(i));

    
    % calculate the difference between top and bottom 50%
    siReal = sort(siReal);

    halfInd = floor(length(siReal)/2);
    
    bottomAvg(i) = mean(siReal(1:halfInd), 'omitnan');
    topAvg(i) = mean(siReal(halfInd+1:end),'omitnan');
    
    difference(i) = topAvg(i) - bottomAvg(i);

    
    
    siShuf = sort(siShuf,1);
    bottomHalfShuf = mean(siShuf(1:halfInd,:),1, 'omitnan');
    topHalfShuf = mean(siShuf(halfInd+1:end,:),1, 'omitnan');
    difShuf = topHalfShuf - bottomHalfShuf;
 
    p(i) = calcP(difference(i), difShuf, 'bigger');
    
    
    sprintf("completed %d of %d",i,length(tVals))
end




%%
clf
hold on

tVals = 1 - tVals;
plot(tVals,topAvg, '-b')
plot(tVals,bottomAvg, '-b')
plot(tVals,difference, '-r')
plot(tVals,p, '--k')
legend('top avg', 'bottom avg', 'difference', 'p-value')


end