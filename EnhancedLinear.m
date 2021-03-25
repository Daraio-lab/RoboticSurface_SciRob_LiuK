function [Rspr, Kspr, Espr] = EnhancedLinear(he,h0,kpi,L0,limlft,limrht)
limlft = limlft/180*pi;
limrht = limrht/180*pi; 
% limlft: theta_1: left partition point
% limrht: theta_2: right partition point
if numel(limlft)==1, limlft = limlft*ones(size(he)); end;
if numel(limrht)==1, limrht = limrht*ones(size(he)); end;
partl = pi./limlft;    partr = pi./(2*pi-limrht);
if numel(kpi)==1, kpi = kpi*ones(size(he)); end;
Rspr = zeros(size(he)); Kspr = Rspr; 

Lind = he<limlft; Rind = he>limrht; Mind = ~(Lind|Rind);
Rspr(Lind) = kpi(Lind).*real(limlft(Lind)-h0(Lind))+kpi(Lind).*tan(partl(Lind)/2.*(he(Lind)-limlft(Lind)))./(partl(Lind)/2);
Kspr(Lind) = kpi(Lind).*sec(partl(Lind)/2.*(he(Lind)-limlft(Lind))).^2;
Rspr(Rind) = kpi(Rind).*real(limrht(Rind)-h0(Rind))+kpi(Rind).*tan(partr(Rind)/2.*(he(Rind)-limrht(Rind)))./(partr(Rind)/2);
Kspr(Rind) = kpi(Rind).*sec(partr(Rind)/2.*(he(Rind)-limrht(Rind))).^2;
Rspr(Mind) = kpi(Mind).*real(he(Mind)-h0(Mind));
Kspr(Mind) = kpi(Mind);
Rspr = L0.*Rspr; Kspr = L0.*Kspr;

if nargout>2
    Espr = zeros(size(he));
    Espr(Lind) = 0.5*kpi(Lind).*real(h0(Lind)-limlft(Lind)).^2+kpi(Lind).*real(h0(Lind)-limlft(Lind)).*(limlft(Lind)-he(Lind))-4*kpi(Lind)./partl(Lind).^2.*log(abs(cos(partl(Lind)/2.*(limlft(Lind)-he(Lind)))));
    Espr(Rind) = 0.5*kpi(Rind).*real(limrht(Rind)-h0(Rind)).^2+kpi(Rind).*real(limrht(Rind)-h0(Rind)).*(he(Rind)-limrht(Rind))-4*kpi(Rind)./partr(Rind).^2.*log(abs(cos(partr(Rind)/2.*(he(Rind)-limrht(Rind)))));
    Espr(Mind) = 0.5*kpi(Mind).*real(he(Mind)-h0(Mind)).^2;
    Espr = L0.*Espr;
end
