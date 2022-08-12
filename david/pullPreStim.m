function [pre, stim] = pullPreStim(g,t)
isPre = t < 0;
isStim = t > .1 & t < 1.1;

pre = g;
pre(~isPre) = nan;

stim = g;
stim(~isStim) = nan;
end