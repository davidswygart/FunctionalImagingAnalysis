function err = errMult(num,den,num_err,den_err)
f = abs(num .* den);
other = sqrt((num ./ num_err).^2 + (den ./ den_err).^2);

err = f .* other;
end