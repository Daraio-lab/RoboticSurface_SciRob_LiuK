function [IF,K] = GlobalK_fast_ver(Ui,Node,truss,angles)
Nn = size(Node,1);
Nodenw(:,1) = Node(:,1)+Ui(1:3:end);
Nodenw(:,2) = Node(:,2)+Ui(2:3:end);
Nodenw(:,3) = Node(:,3)+Ui(3:3:end);

%%
eDofb = kron(truss.Bars,3*ones(1,3))+repmat([-2,-1,0],size(truss.Bars,1),2);
du = [Ui(eDofb(:,1:3))-Ui(eDofb(:,4:6))];
Ex = truss.B*Ui./truss.L+0.5*sum(du.^2,2)./(truss.L.^2);
% Ex = truss.B*Ui./truss.L;
[Sx,Et] = truss.CM(Ex); 
Duelem = [du,-du];
Du = sparse(repmat(1:numel(Et),1,6),eDofb(:),Duelem(:),numel(Et),numel(Ui));
Fx = Sx.*truss.A;
IFb = (sum(bsxfun(@times,truss.B,Fx),1)+sum(bsxfun(@times,Du,Fx./truss.L),1))';
IF = IFb;
if nargout > 1
    Kel = truss.B'*sparse(1:numel(Et),1:numel(Et),Et.*truss.A./truss.L)*truss.B;
    K1 = Du'*sparse(1:numel(Et),1:numel(Et),Et.*truss.A./(truss.L.^2))*truss.B+...
         truss.B'*sparse(1:numel(Et),1:numel(Et),Et.*truss.A./(truss.L.^2))*Du;
    K2 = Du'*sparse(1:numel(Et),1:numel(Et),Et.*truss.A./(truss.L.^3))*Du;
    G = sparse([1:numel(Et),1:numel(Et)],reshape(truss.Bars,[],1),[ones(numel(Et),1);-ones(numel(Et),1)],numel(Et),Nn);
    [ia,ja,sa] = find(G'*bsxfun(@times,G,(Fx./truss.L)));
    ik = bsxfun(@plus, 3*(ia-1), 1:3); jk = bsxfun(@plus, 3*(ja-1), 1:3);
    Kg = sparse(ik,jk,repmat(sa,[1,3]),3*Nn,3*Nn);
    Kg = 0.5*(Kg+Kg');
    Kb = (Kel+K1+K2)+Kg;
    K = Kb;
end

%% Cables
if isfield(truss,'Cables')
    eDofc = kron(truss.Cables,3*ones(1,3))+repmat([-2,-1,0],size(truss.Cables,1),2);
    duc = [Ui(eDofc(:,1:3))-Ui(eDofc(:,4:6))];
    if size(eDofc,1)==1
        duc = duc';
    end
    Exc = truss.Bc*Ui./truss.Lc+0.5*sum(duc.^2,2)./(truss.Lc.^2);
    % Ex = truss.B*Ui./truss.L;
    [Sxc,Etc] = truss.CableCM(Exc,truss.prestrain(:,end)); 
    Duelemc = [duc,-duc];
    Duc = sparse(repmat(1:numel(Etc),1,6),eDofc(:),Duelemc(:),numel(Etc),numel(Ui));
    Fxc = Sxc.*truss.Ac;
    IFc = (sum(bsxfun(@times,truss.Bc,Fxc),1)+sum(bsxfun(@times,Duc,Fxc./truss.Lc),1))';
    IF = IF+IFc;
    if nargout > 1
        Kelc = truss.Bc'*sparse(1:numel(Etc),1:numel(Etc),Etc.*truss.Ac./truss.Lc)*truss.Bc;
        K1c = Duc'*sparse(1:numel(Etc),1:numel(Etc),Etc.*truss.Ac./(truss.Lc.^2))*truss.Bc+...
             truss.Bc'*sparse(1:numel(Etc),1:numel(Etc),Etc.*truss.Ac./(truss.Lc.^2))*Duc;
        K2c = Duc'*sparse(1:numel(Etc),1:numel(Etc),Etc.*truss.Ac./(truss.Lc.^3))*Duc;
        Gc = sparse([1:numel(Etc),1:numel(Etc)],reshape(truss.Cables,[],1),[ones(numel(Etc),1);-ones(numel(Etc),1)],numel(Etc),Nn);
        [iac,jac,sac] = find(Gc'*bsxfun(@times,Gc,(Fxc./truss.Lc)));
        ikc = bsxfun(@plus, 3*(iac-1), 1:3); jkc = bsxfun(@plus, 3*(jac-1), 1:3);
        Kgc = sparse(ikc,jkc,repmat(sac,[1,3]),3*Nn,3*Nn);
        Kc = (Kelc+K1c+K2c)+Kgc;
        K = K+Kc;
    end
end

%% Attaching/Gluing bars
if isfield(truss,'GBars')
    eDofg = kron(truss.GBars,3*ones(1,3))+repmat([-2,-1,0],size(truss.GBars,1),2);
    dug = [Ui(eDofg(:,1:3))-Ui(eDofg(:,4:6))];
    Exg = truss.Bg*Ui./truss.Lg+0.5*sum(dug.^2,2)./(truss.Lg.^2);
    % Ex = truss.B*Ui./truss.L;
    [Sxg,Etg] = truss.CMg(Exg); 
    Duelemg = [dug,-dug];
    Dug = sparse(repmat(1:numel(Etg),1,6),eDofg(:),Duelemg(:),numel(Etg),numel(Ui));
    Fxg = Sxg.*truss.Ag;
    IFg = (sum(bsxfun(@times,truss.Bg,Fxg),1)+sum(bsxfun(@times,Dug,Fxg./truss.Lg),1))';
    IF = IF+IFg;
    if nargout > 1
        Kelg = truss.Bg'*sparse(1:numel(Etg),1:numel(Etg),Etg.*truss.Ag./truss.Lg)*truss.Bg;
        K1g = Dug'*sparse(1:numel(Etg),1:numel(Etg),Etg.*truss.Ag./(truss.Lg.^2))*truss.Bg+...
             truss.Bg'*sparse(1:numel(Etg),1:numel(Etg),Etg.*truss.Ag./(truss.Lg.^2))*Dug;
        K2g = Dug'*sparse(1:numel(Etg),1:numel(Etg),Etg.*truss.Ag./(truss.Lg.^3))*Dug;
        Gg = sparse([1:numel(Etg),1:numel(Etg)],reshape(truss.GBars,[],1),[ones(numel(Etg),1);-ones(numel(Etg),1)],numel(Etg),Nn);
        [iag,jag,sag] = find(Gg'*bsxfun(@times,Gg,(Fxg./truss.Lg)));
        ikg = bsxfun(@plus, 3*(iag-1), 1:3); jkg = bsxfun(@plus, 3*(jag-1), 1:3);
        Kgg = sparse(ikg,jkg,repmat(sag,[1,3]),3*Nn,3*Nn);
        Kbg = (Kelg+K1g+K2g)+Kgg;
        K = K+Kbg;
    end
end

%% Rotational Springs
rotspr = [angles.bend;angles.fold];
eDofd = [kron(rotspr,3*ones(1,3))+repmat([-2,-1,0],size(rotspr,1),4)]';
rkj = [Nodenw(rotspr(:,2),:) - Nodenw(rotspr(:,1),:)]';
rij = [Nodenw(rotspr(:,3),:) - Nodenw(rotspr(:,1),:)]';
rkl = [Nodenw(rotspr(:,2),:) - Nodenw(rotspr(:,4),:)]';
rmj = icross(rij,rkj); rnk = icross(rkj,rkl);
dt_rnkrij = sum(rnk.*rij,1);
sgn = ((abs(dt_rnkrij)>1e-8).*sign(dt_rnkrij)+(abs(dt_rnkrij)<=1e-8)*1);
dt_rmjrnk = sum(rmj.*rnk,1);
rmj2 = sum(rmj.^2,1); norm_rmj = sqrt(rmj2); 
rkj2 = sum(rkj.^2,1); norm_rkj = sqrt(rkj2);
rnk2 = sum(rnk.^2,1); norm_rnk = sqrt(rnk2);
he = sgn.*real(acos(dt_rmjrnk./(norm_rmj.*norm_rnk))); 
he(he<0) = 2*pi+he(he<0); 
Rspr = zeros(size(he));  Kspr = Rspr;
[Rsprb, Ksprb] = angles.CMbend(he(1:size(angles.bend,1))',angles.pb0,angles.Kb,truss.L(1:size(angles.bend,1)));
[Rsprf, Ksprf] = angles.CMfold(he(size(angles.bend,1)+1:end)',angles.pf0,angles.Kf,truss.L(size(angles.bend,1)+1:size(angles.bend,1)+size(angles.fold,1)));
Rspr(1:size(angles.bend,1)) = Rsprb;
Kspr(1:size(angles.bend,1)) = Ksprb;
Rspr(size(angles.bend,1)+1:end) = Rsprf;
Kspr(size(angles.bend,1)+1:end) = Ksprf;

dt_rijrkj = sum(rij.*rkj,1);
dt_rklrkj = sum(rkl.*rkj,1);

di = bsxfun(@times,rmj,norm_rkj./rmj2);
dl = -bsxfun(@times,rnk,norm_rkj./rnk2);
dj = bsxfun(@times,di,(dt_rijrkj./rkj2-1))-bsxfun(@times,dl,dt_rklrkj./rkj2);
dk = -bsxfun(@times,di,dt_rijrkj./rkj2)+bsxfun(@times,dl,(dt_rklrkj./rkj2-1));
Jhe_dense = [dj;dk;di;dl];
Jhe = sparse(eDofd(:),repelem(1:numel(he),12),Jhe_dense(:),numel(Ui),numel(he));
IFbf = sum(bsxfun(@times,Jhe,Rspr),2);
IF = IF+IFbf;

if nargout > 1
    dii = -bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rkj,rmj),[3 1 2]))+...
                         permute(bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rkj,rmj),[3 1 2])),[2 1 3]),...
                  permute(norm_rkj./(rmj2.^2),[3 1 2]));

    dij = -bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(rkj,[3 1 2])),...
                  permute(1./(rmj2.*norm_rkj),[3 1 2]))...
          +bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rkj-rij,rmj),[3 1 2]))+...
                                permute(bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rkj-rij,rmj),[3 1 2])),[2 1 3]),...
                  permute(norm_rkj./(rmj2.^2),[3 1 2]));

    dik =  bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(rkj,[3 1 2])),...
                  permute(1./(rmj2.*norm_rkj),[3 1 2]))...
          +bsxfun(@times,bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rij,rmj),[3 1 2]))+...
                             permute(bsxfun(@times,permute(rmj,[1 3 2]),permute(icross(rij,rmj),[3 1 2])),[2 1 3]),...
                      permute(norm_rkj./(rmj2.^2),[3 1 2]));

    % dil = zeros(3,3,numel(he));

    dll = bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkj,rnk),[3 1 2]))+...
                        permute(bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkj,rnk),[3 1 2])),[2 1 3]),...
                 permute(norm_rkj./(rnk2.^2),[3 1 2]));

    dlk = -bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(rkj,[3 1 2])),...
                  permute(1./(rnk2.*norm_rkj),[3 1 2]))...
          -bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkj-rkl,rnk),[3 1 2]))+...
                         permute(bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkj-rkl,rnk),[3 1 2])),[2 1 3]),...
                  permute(norm_rkj./(rnk2.^2),[3 1 2]));

    dlj =  bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(rkj,[3 1 2])),permute(1./(rnk2.*norm_rkj),[3 1 2]))...
          -bsxfun(@times,bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkl,rnk),[3 1 2]))+...
                         permute(bsxfun(@times,permute(rnk,[1 3 2]),permute(icross(rkl,rnk),[3 1 2])),[2 1 3]),...
                  permute(norm_rkj./(rnk2.^2),[3 1 2]));

    dT1jj = bsxfun(@times,bsxfun(@times,rkj,(-1+2*dt_rijrkj./rkj2))-rij,1./(rkj2));
    dT2jj = bsxfun(@times,bsxfun(@times,rkj,(2*dt_rklrkj./rkj2))-rkl,1./(rkj2));
    djj = bsxfun(@times,permute(di,[1 3 2]),permute(dT1jj,[3 1 2]))+...
          bsxfun(@times,dij,permute((dt_rijrkj./rkj2)-1,[3 1 2]))-...
          bsxfun(@times,permute(dl,[1 3 2]),permute(dT2jj,[3 1 2]))-...
          bsxfun(@times,dlj,permute(dt_rklrkj./rkj2,[3 1 2]));

    dT1jk = bsxfun(@times,bsxfun(@times,rkj,(-2*dt_rijrkj./rkj2))+rij,1./(rkj2));
    dT2jk = bsxfun(@times,bsxfun(@times,rkj,(1-2*dt_rklrkj./rkj2))+rkl,1./(rkj2));
    djk = bsxfun(@times,permute(di,[1 3 2]),permute(dT1jk,[3 1 2]))+...
          bsxfun(@times,dik,permute((dt_rijrkj./rkj2)-1,[3 1 2]))-...
          bsxfun(@times,permute(dl,[1 3 2]),permute(dT2jk,[3 1 2]))-...
          bsxfun(@times,dlk,permute(dt_rklrkj./rkj2,[3 1 2]));

    dT1kk = dT2jk; dT2kk = dT1jk;
    dkk = bsxfun(@times,permute(dl,[1 3 2]),permute(dT1kk,[3 1 2]))+...
          bsxfun(@times,dlk,permute(dt_rklrkj./rkj2-1,[3 1 2]))-...
          bsxfun(@times,permute(di,[1 3 2]),permute(dT2kk,[3 1 2]))-...
          bsxfun(@times,dik,permute(dt_rijrkj./rkj2,[3 1 2]));

    Hp = zeros(12,12,numel(he));
    Hp(1:3,1:3,:) = djj; Hp(7:9,7:9,:) = dii; Hp(4:6,4:6,:) = dkk; Hp(10:12,10:12,:) = dll;
    Hp(1:3,4:6,:) = djk;   Hp(4:6,1:3,:) = permute(djk,[2 1 3]);
    Hp(7:9,1:3,:) = dij; Hp(1:3,7:9,:) = permute(dij,[2 1 3]);
    Hp(10:12,1:3,:) = dlj; Hp(1:3,10:12,:) = permute(dlj,[2 1 3]);
    Hp(7:9,4:6,:) = dik;  Hp(4:6,7:9,:) = permute(dik,[2 1 3]);
    Hp(10:12,4:6,:) = dlk;  Hp(4:6,10:12,:) = permute(dlk,[2 1 3]);

    Khe_dense = bsxfun(@times,bsxfun(@times,permute(Jhe_dense,[1 3 2]),permute(Jhe_dense,[3 1 2])),permute(Kspr,[3 1 2]))+...
                bsxfun(@times,Hp,permute(Rspr,[3 1 2]));

    dof_ind1 = permute(repmat(eDofd,[1 1 12]),[1 3 2]);
    dof_ind2 = permute(dof_ind1,[2 1 3]);
    Kbf = sparse(dof_ind1(:),dof_ind2(:),Khe_dense(:),3*Nn,3*Nn);

    K = K+Kbf;
    K = (K+K')/2;
    
    K = K+sparse(1:1:3*Nn,1:1:3*Nn,1e-5*ones(1,3*Nn),3*Nn,3*Nn);
end
end

%--------------------------------------------------------------------------
function c = icross(a,b)
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:)
     a(3,:).*b(1,:)-a(1,:).*b(3,:)
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
end


