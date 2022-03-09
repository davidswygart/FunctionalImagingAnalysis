function [SI, err] = calcSI_DFoF(A,a,B,b,sA,sa,sB,sb)
% calcSI = @(A,a,B,b)(A.*b - B.*a) ./ (A.*b + B.*a - 2*b.*a); %SI using dF/F
% num,den,num_err,den_err
Ab = A .* b;
Ba = B .* a;
ab2 = a .* b .* 2;

Ab_err = errMult(A,b,sA,sb);
Ba_err = errMult(B,a,sB,sa);
ab2_err = errMult(a,b,sa,sb) * 2;

Ab_Ba = Ab - Ba;
Ab_Ba_err = sqrt(Ab_err.^2 + Ba_err.^2);

Ab_Ba_2ab = Ab + Ba - ab2;
Ab_Ba_2ab_err = sqrt(Ab_err.^2 + Ba_err.^2 + ab2_err.^2);

SI = Ab_Ba ./ Ab_Ba_2ab;
err = errDiv(Ab_Ba, Ab_Ba_2ab, Ab_Ba_err, Ab_Ba_2ab_err);
end