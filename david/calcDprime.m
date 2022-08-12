function [dprime, difs] = calcDprime(si,err)
si = si(:);
err = err(:);

d = abs(si - si');
avgErr = (err + err') / 2;

dprime = d ./ avgErr;

isU = logical(tril(ones(size(dprime)),-1));

difs = d(isU);
dprime = dprime(isU);
end