function [Rspr, Kspr, Espr] = SuperLinearBend(he,h0,kp,L0)
Rspr = sign(he-h0).*abs(he-h0).^(4/3).*kp.*L0;
Kspr = max(4/3*abs(he-h0).^(1/3).*kp,0.001*kp).*L0;
if nargout>2
    Espr = 3/7*abs(he-h0).^(7/3).*kp.*L0;
end