function [medSI,allSI] = calcAndPlotSi(small,big, varargin)
%% calculate si
allSI = big ./ small;
medSI = median(allSI,3, 'omitnan');


if ~isempty(varargin) % plot SI
    clf
    imagesc(medSI)
    c = colorbar;
    c.Label.String = 'Big / Small';

    caxis([-.5,1.5])
    axis image
    title('median suppression')
end

%madSI = mad(allSI,1,3);
% 
% %% Show error of SI (MAD)
% figure(8)
% clf
% imagesc(madSI)
% c = colorbar;
% c.Label.String = 'median absolute deviation';
% 
% caxis([-.5,1.5])
% axis image
% title('si error (MAD)')
end