set(groot,'defaultLineLineWidth',2.0)
%% get response and null values
smallDf = eResp.dF(:,:,smallInds);
bigDf = eResp.dF(:,:,bigInds);

nullDf1 = eResp.dF(:,:,nullInds1);
nullDf2 = eResp.dF(:,:,nullInds2);

%% plot si
figure(14)
clf
[medSI,allSI] = calcAndPlotSi(smallDf,bigDf, 1);
saveas(gcf,[image,'_7_si_real.png'])

%% plot null si
figure(15)
clf
[siN,siNall] = calcAndPlotSi(nullDf1,nullDf2, 1);
saveas(gcf,[image,'_8_si_null.png'])
% % plot null si error
% figure(16)
% clf
siN_mad = mad(nullDf2 ./ nullDf1, 1,3);
% imagesc(siN_mad)
% colorbar
% caxis([0,2])
%% correlations to Null data
correlations = [
corNan(siN,smallSNR)
corNan(siN,smallDfof)
corNan(siN,median(smallDf,3))
corNan(siN,rawG)
corNan(siN,siN_mad)];

X = categorical({'SNR','dF/F','dF','raw','si MAD'});
X = reordercats(X,{'SNR','dF/F','dF','raw','si MAD'});


%% Null 1 p values (using smallSNR)
figure(8)
[p,v,tVals,pixPerc] = pByThresh(nullDf1, nullDf2, siN);
title('Null 1: using si null')
%saveas(gcf,[image,'_10_Null_SNR.png'])
%% Null of real values (using smallSNR)
figure(9)
[p,v,tVals,pixPerc] = pByThresh(smallDf, bigDf, siN);
title('real: using si null')
%saveas(gcf,[image,'_11_Real_SNR.png'])


% figure(7)
% clf
% bar(X,correlations)
% ylabel('correlation')
% saveas(gcf,[image,'_9_Null_Correlation.png'])
% 
% 
% %% Null 1 p values (using smallSNR)
% figure(8)
% [p,v,tVals,pixPerc] = pByThresh(nullDf1, nullDf2, smallSNR);
% title('Null 1: using SNR')
% saveas(gcf,[image,'_10_Null_SNR.png'])
% %% Null of real values (using smallSNR)
% figure(9)
% [p,v,tVals,pixPerc] = pByThresh(smallDf, bigDf, smallSNR);
% title('real: using SNR')
% saveas(gcf,[image,'_11_Real_SNR.png'])
% 
% 
% %% Null 1 p values (using smallDfof)
% figure(10)
% [p,v,tVals,pixPerc] = pByThresh(nullDf1, nullDf2, smallDfof);
% title('Null 1: using df / f')
% saveas(gcf,[image,'_12_Null_dFoF.png'])
% %% Null of real values (using smallDfof)
% figure(11)
% [p,v,tVals,pixPerc] = pByThresh(smallDf, bigDf, smallDfof);
% title('real: using df / f')
% saveas(gcf,[image,'_13_Real_dFoF.png'])
% 
% 
% %% Null 1 p values (using raw)
% figure(12)
% [p,v,tVals,pixPerc] = pByThresh(nullDf1, nullDf2, rawG);
% title('Null 1: using SNR')
% saveas(gcf,[image,'_14_Null_Raw.png'])
% %% Null of real values (using raw)
% figure(13)
% [p,v,tVals,pixPerc] = pByThresh(smallDf, bigDf, rawG);
% title('real: using SNR')
% saveas(gcf,[image,'_15_Real_Raw.png'])
% 
% 
