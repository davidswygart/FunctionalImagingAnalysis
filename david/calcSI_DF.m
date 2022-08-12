function [SI, err] = calcSI_DF(A,a,B,b,sA,sa,sB,sb)
%calcSI = @(A,a,B,b)((A-a-B+b) ./ (A-a+B-b)); %SI using dF
% actual numerator and denomenator values
num = A-a-B+b;
den = A-a+B-b;

% error for numerator and denomenator values
num_err = sqrt(sA.^2 + sa.^2 + sB.^2 + sb.^2);
den_err = num_err;

SI = num ./ den;
err = errDiv(num,den,num_err,den_err);
end