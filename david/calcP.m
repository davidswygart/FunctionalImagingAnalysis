function p = calcP(real,shuf,direction)

shuf(isnan(shuf)) = [];

if strcmp(direction,'smaller')
    p = mean(real > shuf);
elseif strcmp(direction,'bigger')
    p = mean(real < shuf);
end

% clf
% histogram(shuf)
% hold on
% yLims = ylim;
% plot([real,real],yLims, '--k')
% text(real,yLims(2)*.8,['  p = ',num2str(p)])

end