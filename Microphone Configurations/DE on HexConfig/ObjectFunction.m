function z = ObjectFunction(array, n_microphones, alpha_min, alpha_max)

% Eq. (52) of NLR-TP-2010-549

myEq = 0;

for i = 1:(n_microphones-1)
   a = array(i)-array(i+1:n_microphones);
   b = array(i+n_microphones)-array(i+n_microphones+1:2*n_microphones);
   dis = sqrt(a.^2+b.^2);
   myEq = myEq + sum(2*(alpha_max*-1*besselj(1, alpha_max.*dis) - ...
            alpha_min*-1*besselj(1, alpha_min.*dis)) ./ dis);
end

z = (pi/n_microphones^2) * (n_microphones*(alpha_max^2-alpha_min^2) - 2*myEq);