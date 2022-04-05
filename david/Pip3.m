set(groot,'defaultLineLineWidth',2.0)
%% get response and null values
smallDf = eResp.dF(:,:,smallInds);
bigDf = eResp.dF(:,:,bigInds);

nullDf1 = eResp.dF(:,:,nullInds1);
nullDf2 = eResp.dF(:,:,nullInds2);

% %% plot si
% figure(14)
% clf
% [medSI,allSI] = calcAndPlotSi(smallDf,bigDf, 1);
% saveas(gcf,[image,'_7_si_real.png'])

%% plot null si
figure(15)
clf
[siN,siNall] = calcAndPlotSi(nullDf1,nullDf2, 1);
caxis([.6 1.1])
title(['si: ', sizeStr, '/ ', sizeStr])
saveas(gcf,[image,'_7_si_null.png'])


% %% 
% figure(16)
% clf
% diffByThresh(nullDf1,nullDf2, siN);
% 
% %% 
% figure(17)
% clf
% diffByThresh(smallDf, bigDf, siN);
% 





