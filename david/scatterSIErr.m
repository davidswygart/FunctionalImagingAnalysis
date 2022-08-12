function scatterSIErr(si,err)
si = si(:);
err = err(:);

figure(1)
clf

scatterhist(si,err);
xlabel('SI')
ylabel('std')



end