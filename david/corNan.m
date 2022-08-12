function corrVal = corNan(a,b)
a = a(:);
b = b(:);

isbad = isnan(a) | isnan(b) | a==inf | a==-inf | b==-inf | b==-inf;

a(isbad) = [];
b(isbad) = [];

corrVal = corr(a,b);
end