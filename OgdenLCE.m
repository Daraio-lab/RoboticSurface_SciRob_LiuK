function [Sx, Ct, Wb] = OgdenLCE(Ex, Pstr0)
% Ogden hyperelastic constitutive model for bar elements
p = [0.4183, 6.1965, 0.8086, 0.8382, -0.4900, 5.9338];
% alfa = [3,1]; % Linear
pstr = real(sqrt(2*Ex+1))./Pstr0;

Ct = p(1)*((p(2)-2)*pstr.^(p(2)-4)+(0.5*p(2)+2)*pstr.^(-0.5*p(2)-4))+...
     p(3)*((p(4)-2)*pstr.^(p(4)-4)+(0.5*p(4)+2)*pstr.^(-0.5*p(4)-4))+...
     p(5)*((p(6)-2)*pstr.^(p(6)-4)+(0.5*p(6)+2)*pstr.^(-0.5*p(6)-4));
 
Ct(pstr<1) = 0;
 
Sx = p(1)*(pstr.^(p(2)-2)-pstr.^(-0.5*p(2)-2))+...
     p(3)*(pstr.^(p(4)-2)-pstr.^(-0.5*p(4)-2))+...
     p(5)*(pstr.^(p(6)-2)-pstr.^(-0.5*p(6)-2));

Sx(pstr<1) = 0;
 
if nargout>2
    Wb = p(1)/p(2)*(pstr.^p(2)+2*pstr.^(-0.5*p(2))-3)+...
         p(3)/p(4)*(pstr.^p(4)+2*pstr.^(-0.5*p(4))-3)+...
         p(5)/p(6)*(pstr.^p(6)+2*pstr.^(-0.5*p(6))-3);
    Wb(pstr<1) = 0;
end
    