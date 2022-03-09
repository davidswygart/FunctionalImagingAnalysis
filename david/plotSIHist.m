function plotSIHist(dF_SI, dFoF_SI, center)
    edges = linspace(center-1,center+1,10);
    edges(1) = -inf;
    edges(end) = inf;

    figure(1)
    clf
    subplot(2,1,1)
    v = dF_SI;
    histogram(dF_SI,edges)
    title('using dF')
    ylabel(['std = ' num2str(std(v))])
    xlabel(['within 1 = ', num2str(mean(abs(v-center) < .5))])

    subplot(2,1,2)
    v = dFoF_SI;
    histogram(v,edges)
    title('using dF/F')
    ylabel(['std = ' num2str(std(v))])
    xlabel(['within 1 = ', num2str(mean(abs(v-center) < .5))])
end