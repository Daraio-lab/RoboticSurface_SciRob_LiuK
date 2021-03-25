function STAT = PostProcess(U_his,truss,angles)
%% Get Data
Exbar = zeros(size(truss.Bars,1),size(U_his,2));
if isfield(truss,'Cables')
    Excable = zeros(size(truss.Cables,1),size(U_his,2)); 
    Sx_cable = zeros(size(truss.Cables,1),size(U_his,2));
    Wc = zeros(size(truss.Cables,1),size(U_his,2));
end
if isfield(truss,'GBars')
    Exgbar = zeros(size(truss.GBars,1),size(U_his,2));
end
FdAngle = zeros(size(angles.fold,1),size(U_his,2)); 
BdAngle = zeros(size(angles.bend,1),size(U_his,2));
for icrm=1:size(U_his,2)
    Ui = U_his(:,icrm);
    Nodenw = truss.Node;
    Nodenw(:,1) = truss.Node(:,1)+Ui(1:3:end);
    Nodenw(:,2) = truss.Node(:,2)+Ui(2:3:end);
    Nodenw(:,3) = truss.Node(:,3)+Ui(3:3:end);
    
    eDofb = kron(truss.Bars,3*ones(1,3))+repmat([-2,-1,0],size(truss.Bars,1),2);
    du = [Ui(eDofb(:,1:3))-Ui(eDofb(:,4:6))];
    Exbar(:,icrm) = truss.B*Ui./truss.L+0.5*sum(du.^2,2)./(truss.L.^2);
    
    if isfield(truss,'Cables')
        eDofc = kron(truss.Cables,3*ones(1,3))+repmat([-2,-1,0],size(truss.Cables,1),2);
        duc = [Ui(eDofc(:,1:3))-Ui(eDofc(:,4:6))];
        if size(eDofc,1)==1
            duc = duc';
        end
        Excablei = truss.Bc*Ui./truss.Lc+0.5*sum(duc.^2,2)./(truss.Lc.^2);
        [Sx_cablei, ~, Wci] = truss.CableCM(Excablei, truss.prestrain(icrm+1));
        Excable(:,icrm) = Excablei;
        Sx_cable(:,icrm) = Sx_cablei;
        Wc(:,icrm) = Wci;
    end
    
    if isfield(truss,'GBars')
        eDofg = kron(truss.GBars,3*ones(1,3))+repmat([-2,-1,0],size(truss.GBars,1),2);
        dug = [Ui(eDofg(:,1:3))-Ui(eDofg(:,4:6))];
        Exgbar(:,icrm) = truss.Bg*Ui./truss.Lg+0.5*sum(dug.^2,2)./(truss.Lg.^2);
    end

    for del = 1:size(angles.bend,1)
        bend = angles.bend(del,:);
        BdAngle(del,icrm) = FoldKe(Nodenw,bend);
    end

    for fel = 1:size(angles.fold,1)
        fold = angles.fold(fel,:);
        FdAngle(fel,icrm) = FoldKe(Nodenw,fold);
    end
end

%% Interpret Data
[Sx_bar, ~, Wb] = truss.CM(Exbar);
if isfield(truss,'GBars')
    [Sx_gbar, ~, Wg] = truss.CMg(Exgbar);
end
Rspr_fd = zeros(size(FdAngle)); Efold = Rspr_fd;
Rspr_bd = zeros(size(BdAngle)); Ebend = Rspr_bd;
for i = 1:size(U_his,2)
    [Rspr_fdi, ~, Efoldi] = angles.CMfold(FdAngle(:,i),angles.pf0,angles.Kf,truss.L((size(angles.bend,1)+1):(size(angles.bend,1)+size(angles.fold,1))));
    [Rspr_bdi, ~, Ebendi] = angles.CMbend(BdAngle(:,i),angles.pb0,angles.Kb,truss.L(1:size(angles.bend,1)));
    Rspr_fd(:,i) = Rspr_fdi; Efold(:,i) = Efoldi;
    Rspr_bd(:,i) = Rspr_bdi; Ebend(:,i) = Ebendi;
end

STAT.bar.Ex = Exbar; 
STAT.bar.Sx = Sx_bar;    
STAT.bar.USi = diag(truss.L.*truss.A)*Wb;
STAT.bar.US = sum(STAT.bar.USi,1);

if isfield(truss,'Cables')
    STAT.cable.Ex = Excable; 
    STAT.cable.Sx = Sx_cable;    
    STAT.cable.USi = diag(truss.Lc.*truss.Ac)*Wc;
    STAT.cable.US = sum(STAT.cable.USi,1);
end

if isfield(truss,'GBars')
    STAT.gbar.Ex = Exgbar; 
    STAT.gbar.Sx = Sx_gbar;    
    STAT.gbar.USi = diag(truss.Lg.*truss.Ag)*Wg;
    STAT.gbar.US = sum(STAT.gbar.USi,1);
end

STAT.fold.Angle = FdAngle;
STAT.fold.RM = Rspr_fd;
STAT.fold.UFi = Efold;
STAT.fold.UF = sum(STAT.fold.UFi,1);

STAT.bend.Angle = BdAngle;
STAT.bend.RM = Rspr_bd;
STAT.bend.UBi = Ebend;
STAT.bend.UB = sum(STAT.bend.UBi,1);

STAT.PE = STAT.bar.US+STAT.fold.UF+STAT.bend.UB;

if isfield(truss,'Cables')
    STAT.PE = STAT.PE+STAT.cable.US;
end
if isfield(truss,'GBars')
    STAT.PE = STAT.PE+STAT.gbar.US;
end
