function [truss, angles, AnalyInputOpt] = PrepareData(Node,Panel,Supp,Load,AnalyInputOpt)
% AnalyInputOpt contains the following fields:
%     'ModelType': 'N4B5' or 'N5B8'
%     'MaterCalib': 'auto' or 'manual'
%     'BarCM': Constitutive model of bar elements
%     'Abar': Area of bars (manual mode)
%     'RotSprBend': Constitutive model of bending elements (manual mode)
%     'RotSprFold': Constitutive model of folding elements (manual mode)
%     'Kb': Initial modulus of bending hinges
%     'Kf': Initial modulus of folding hinges
%     'ModElastic': Modulus of Elasticy of bar material
%     'Poisson': Poisson's ratio
%     'Thickness': Panel thickness
%     'LScaleFactor': Ratio of length scale factor (Ls) over hinge length LF
%     'ZeroBend': 'AsIs' - as the geometry is defined, 
%                 'Flat' - enforce neutral angle at pi (flat configuration)
%                  value - provide specific values, scalar for uniform
%                  assignment, vector for differential assignment 
%     'LoadType': 'Force' or 'Displacement'
%     'Load': Loading condition
%     'AdaptiveLoad': Function handle for adaptive loading
%     'InitialLoadFactor': Initial load factor for MGDCM ('Force' mode)
%     'MaxIcr': Maximum number of increments ('Force' mode)
%     'DispStep': Number of displacement increments ('Displacement' mode)
%     'StopCriterion': Function handle for specific stopping criterion.
%     'Cables': auxiliary 1D elements
%     'CableCM': Constitutive model of auxiliary 1D elements
%     'CableA': areas of auxiliary 1D elements
%     'GBars': auxiliary 1D elements
%     'GBarCM': Constitutive model of auxiliary 1D elements
%     'GBarA': areas of auxiliary 1D elements

if ~isfield(AnalyInputOpt,'StopCriterion')
    AnalyInputOpt.StopCriterion = @(Node,U,icrm)false;
end

if ~isfield(AnalyInputOpt,'ModelType')
    AnalyInputOpt.ModelType = 'N5B8';
end

if ~isfield(AnalyInputOpt,'LoadType')
    AnalyInputOpt.LoadType = 'Force';
end

if ~isfield(AnalyInputOpt,'MaxIcr')
    AnalyInputOpt.MaxIcr = 100;
end

if ~isfield(AnalyInputOpt,'InitialLoadFactor')
    AnalyInputOpt.InitialLoadFactor = 0.01;
end

if strcmpi(AnalyInputOpt.ModelType,'N4B5')
    % N4B5 model only accepts manually input for material properties
    [Bend, Node] = findbend(Panel, Node, AnalyInputOpt.ModelType);
    % Find bending hinges
    [Fold, Bdry, Trigl] = findfdbd(Panel,Bend);
    % Find folding hinges and boundaries, return final triangulation
    Bars = [Bend(:,1:2);Fold(:,1:2);Bdry];
    % Define bar elements
    [B, L] = dirc3d(Node,Bars);
    if isfield(AnalyInputOpt,'Cables')
        [Bc, Lc] = dirc3d(Node,AnalyInputOpt.Cables);
    end
    if isfield(AnalyInputOpt,'GBars')
        [Bg, Lg] = dirc3d(Node,AnalyInputOpt.GBars);
    end
    if size(Supp,1) == 0
        rs = []; 
    else
        rs = [reshape([Supp(:,1)*3-2,Supp(:,1)*3-1,Supp(:,1)*3]',[],1),...
              reshape(Supp(:,2:4)',[],1)];
        rs(rs(:,2)==0,:)=[]; rs = rs(:,1);
    end 
    
    Abar = getfieldvalues(AnalyInputOpt,'Abar',0.1);
    if numel(Abar)==1
        Abar = Abar*ones(size(Bars,1),1);
    end
    
    if isfield(AnalyInputOpt,'Cables')
        ACable = getfieldvalues(AnalyInputOpt,'CableA',0.1);
        if numel(ACable)==1
            ACable = ACable*ones(size(AnalyInputOpt.Cables,1),1);
        end
    end
    
    if isfield(AnalyInputOpt,'GBars')
        AGBar = getfieldvalues(AnalyInputOpt,'GBarA',0.1);
        if numel(AGBar)==1
            AGBar = AGBar*ones(size(AnalyInputOpt.GBars,1),1);
        end
    end

    pf0 = zeros(size(Fold,1),1); 
    for i = 1:size(Fold,1), pf0(i) = FoldKe(Node,Fold(i,:)); end

    ZeroBend = getfieldvalues(AnalyInputOpt,'ZeroBend','AsIs'); 
    pb0 = zeros(size(Bend,1),1); 
    if isnumeric(ZeroBend)
        pb0 = pb0+ZeroBend;
    elseif strcmpi(ZeroBend,'Flat')
        pb0 = pb0+pi;    
    elseif strcmpi(ZeroBend,'AsIs')
        for i = 1:size(Bend,1), pb0(i) = FoldKe(Node,Bend(i,:)); end
    end
    
    kpb = getfieldvalues(AnalyInputOpt,'Kb',0.1);
    if (size(kpb,1)==1)&&(size(Bend,1)>1)
        kpb = repmat(kpb,size(Bend,1),1);
    end
    
    kpf = getfieldvalues(AnalyInputOpt,'Kf',0.1*kpb(1));
    if (size(kpf,1)==1)&&(size(Fold,1)>1)
        kpf = repmat(kpf,size(Fold,1),1);
    end
    
    if ~isempty(Load)
        m = size(Node,1);
        FD = zeros(3*m,1);
        indp = Load(:,1);
        FD(3*indp-2) = Load(:,2); 
        FD(3*indp-1) = Load(:,3); 
        FD(3*indp) = Load(:,4);
        AnalyInputOpt.Load = FD;
    end

    truss.CM = getfieldvalues(AnalyInputOpt,'BarCM',@(Ex)Ogden(Ex, 1e4));
    if isfield(AnalyInputOpt,'Cables')
        truss.CableCM = getfieldvalues(AnalyInputOpt,'CableCM',@(Ex)Ogden(Ex, 12));
        truss.Cables = AnalyInputOpt.Cables;
        truss.Bc = Bc; 
        truss.Lc = Lc;
        truss.Ac = ACable;
        truss.prestrain = ones(size(ACable,1),1);
    end
    if isfield(AnalyInputOpt,'GBars')
        truss.GBarCM = getfieldvalues(AnalyInputOpt,'GBarCM',@(Ex)Ogden(Ex, 120));
        truss.GBars = AnalyInputOpt.GBars;
        truss.Bg = Bg; 
        truss.Lg = Lg;
        truss.Ag = AGBar;
    end
    truss.Node = Node;
    truss.Bars = Bars;
    truss.Trigl = Trigl;
    truss.B = B; 
    truss.L = L;
    truss.FixedDofs = unique(rs);
    truss.A = Abar;
    angles.CMbend = getfieldvalues(AnalyInputOpt,'RotSprBend',@(he,h0,kb,L0)EnhancedLinear(he,h0,kb,L0,45,315));
    angles.CMfold = getfieldvalues(AnalyInputOpt,'RotSprFold',@(he,h0,kf,L0)EnhancedLinear(he,h0,kf,L0,45,315));
    angles.fold = Fold;
    angles.bend = Bend;
    angles.Kb = kpb;
    angles.Kf = kpf;
    angles.pf0 = pf0;
    angles.pb0 = pb0;
    angles.Panel = Panel;
    
elseif strcmpi(AnalyInputOpt.ModelType,'N5B8')
    [Bend, Node, panelctr] = findbend(Panel, Node, AnalyInputOpt.ModelType);
    [Fold, Bdry, Trigl] = findfdbd(Panel,Bend);
    Bars = [Bend(:,1:2);Fold(:,1:2);Bdry];
    [B, L] = dirc3d(Node,Bars);
    if isfield(AnalyInputOpt,'Cables')
        [Bc, Lc] = dirc3d(Node,AnalyInputOpt.Cables);
    end
    if isfield(AnalyInputOpt,'GBars')
        [Bg, Lg] = dirc3d(Node,AnalyInputOpt.GBars);
    end
    if size(Supp,1) == 0
        rs = []; 
    else
        rs = [reshape([Supp(:,1)*3-2,Supp(:,1)*3-1,Supp(:,1)*3]',[],1),...
              reshape(Supp(:,2:4)',[],1)];
        rs(rs(:,2)==0,:)=[]; rs = rs(:,1);
    end
    
    pf0 = zeros(size(Fold,1),1); 
    for i = 1:size(Fold,1), pf0(i) = FoldKe(Node,Fold(i,:)); end

    ZeroBend = getfieldvalues(AnalyInputOpt,'ZeroBend','AsIs'); 
    pb0 = zeros(size(Bend,1),1); 
    if isnumeric(ZeroBend)
        pb0 = pb0+ZeroBend;
    elseif strcmpi(ZeroBend,'Flat')
        pb0 = pb0+pi;    
    elseif strcmpi(ZeroBend,'AsIs')
        for i = 1:size(Bend,1), pb0(i) = FoldKe(Node,Bend(i,:)); end
    end
    
    if ~isempty(Load)
        m = size(Node,1);
        FD = zeros(3*m,1);
        indp = Load(:,1);
        FD(3*indp-2) = Load(:,2); 
        FD(3*indp-1) = Load(:,3); 
        FD(3*indp) = Load(:,4);
        AnalyInputOpt.Load = FD;
    end
    
    if isfield(AnalyInputOpt,'Cables')
        ACable = getfieldvalues(AnalyInputOpt,'CableA',0.1);
        if numel(ACable)==1
            ACable = ACable*ones(size(AnalyInputOpt.Cables,1),1);
        end
    end
    
    if isfield(AnalyInputOpt,'GBars')
        AGBar = getfieldvalues(AnalyInputOpt,'GBarA',0.1);
        if numel(AGBar)==1
            AGBar = AGBar*ones(size(AnalyInputOpt.GBars,1),1);
        end
    end
    
    if isfield(AnalyInputOpt,'Cables')
        truss.CableCM = getfieldvalues(AnalyInputOpt,'CableCM',@(Ex)Ogden(Ex, 12));
        truss.Cables = AnalyInputOpt.Cables;
        truss.Bc = Bc; 
        truss.Lc = Lc;
        truss.Ac = ACable;
        truss.prestrain = ones(size(ACable,1),1);
    end
    if isfield(AnalyInputOpt,'GBars')
        truss.CMg = getfieldvalues(AnalyInputOpt,'GBarCM',@(Ex)Ogden(Ex, 120));
        truss.GBars = AnalyInputOpt.GBars;
        truss.Bg = Bg; 
        truss.Lg = Lg;
        truss.Ag = AGBar;
    end
    truss.Node = Node;
    truss.Bars = Bars;
    truss.Trigl = Trigl;
    truss.B = B; 
    truss.L = L;
    truss.FixedDofs = unique(rs);
    angles.fold = Fold;
    angles.bend = Bend;
    angles.pf0 = pf0;
    angles.pb0 = pb0;
    angles.Panel = Panel;
    MaterCalib = getfieldvalues(AnalyInputOpt,'MaterCalib','auto');
    if strcmpi(MaterCalib,'manual')
        truss.CM = getfieldvalues(AnalyInputOpt,'BarCM',@(Ex)Ogden(Ex, 1e4));
        
        Abar = getfieldvalues(AnalyInputOpt,'Abar',0.1);
        if numel(Abar)==1
            Abar = Abar*ones(size(Bars,1),1);
        end

        kpb = getfieldvalues(AnalyInputOpt,'Kb',0.1);
        if (size(kpb,1)==1)&&(size(Bend,1)>1)
            kpb = repmat(kpb,size(Bend,1),1);
        end

        kpf = getfieldvalues(AnalyInputOpt,'Kf',0.1*kpb(1));
        if (size(kpf,1)==1)&&(size(Fold,1)>1)
            kpf = repmat(kpf,size(Fold,1),1);
        end
        
        truss.A = Abar;
        angles.CMbend = getfieldvalues(AnalyInputOpt,'RotSprBend',@(he,h0,kb,L0)EnhancedLinear(he,h0,kb,L0,45,315));
        angles.CMfold = getfieldvalues(AnalyInputOpt,'RotSprFold',@(he,h0,kf,L0)EnhancedLinear(he,h0,kf,L0,45,315));
        angles.Kb = kpb;
        angles.Kf = kpf;
    else
        EY = getfieldvalues(AnalyInputOpt,'ModElastic',1e9);
        nv = getfieldvalues(AnalyInputOpt,'Poisson',0.33);
        thck = getfieldvalues(AnalyInputOpt,'Thickness',0.15e-2);
        truss.CM = getfieldvalues(AnalyInputOpt,'BarCM',@(Ex)Ogden(Ex, EY));
        Abar = zeros(size(Bars,1),1);
        kpbt = zeros(size(Bend,1),1);
        kpft = zeros(size(Fold,1),1);
        G = EY*thck^3/(12*(1-nv^2));
        Lf = L((numel(kpbt)+1):(numel(kpbt)+numel(kpft)));
        Ls = getfieldvalues(AnalyInputOpt,'LScaleFactor',2*sum(Lf)/numel(Lf));
        Kl = 1/Ls*G;  Km = 0.55*G*(Lf/thck).^(1/3);
        kpft = (1./(1./Kl+1./Km))./Lf;

        for j=1:length(panelctr)       
                [Abarj, kpbj] = GetMaterial(Node,Panel{j},panelctr(j),EY,nv,thck,Fold,Bend,Bdry);
                Abar(Abarj(:,1)) = Abar(Abarj(:,1))+Abarj(:,2);
                if ~isempty(kpbj)
                    kpbt(kpbj(:,1)) = kpbt(kpbj(:,1))+kpbj(:,2);
                end
        end
            
        truss.A = Abar;
        angles.CMbend = getfieldvalues(AnalyInputOpt,'RotSprBend',@(he,h0,kb,L0)EnhancedLinear(he,h0,kb,L0,30,330)); 
        angles.CMfold = getfieldvalues(AnalyInputOpt,'RotSprFold',@(he,h0,kf,L0)EnhancedLinear(he,h0,kf,L0,30,330));
        
        kpb = getfieldvalues(AnalyInputOpt,'Kb',kpbt);
        if (size(kpb,1)==1)&&(size(Bend,1)>1)
            kpb = repmat(kpb,size(Bend,1),1);
        end

        kpf = getfieldvalues(AnalyInputOpt,'Kf',kpft);
        if (size(kpf,1)==1)&&(size(Fold,1)>1)
            kpf = repmat(kpf,size(Fold,1),1);
        end
        angles.Kb = kpb;
        angles.Kf = kpf;
    end
else
    error('Model type not supported!')
end
end


%--------------------------------------------------------------------------
function [bend, Node, panelctr] = findbend(Panel, Node, ModelType)
if ~strcmp(ModelType,'N5B8')   %'N4B5'
    bend = nan(3*length(Panel),4);
    cntb = 0;
    for i = 1:length(Panel)
        if numel(Panel{i}) == 4
            L1 = norm((Node(Panel{i}(1),:)-Node(Panel{i}(3),:)));
            L2 = norm((Node(Panel{i}(4),:)-Node(Panel{i}(2),:)));
            if L1>L2, lclbend = [2,4,1,3];
            else lclbend = [1,3,2,4]; 
            end
            cntb = cntb+1;
            bend(cntb,:) = Panel{i}(lclbend);
        elseif numel(Panel{i}) > 4
            BiEdge = DividePolygon(Node(Panel{i},:));
            if (numel(Panel{i})-3)==size(BiEdge,1)
                bendlcl = nan(size(BiEdge,1),4);
            else
                error('Impossible Triangluation!')
            end
            pp = [1:numel(Panel{i})]';
            TotalLinks = [[circshift(pp,1,1),pp];...
                          BiEdge];
            Comm = sparse(numel(pp),size(TotalLinks,1));
            for j=1:size(TotalLinks,1), Comm(TotalLinks(j,:),j) = true; end
            Gp = Comm*Comm'; 
            for k = 1:size(BiEdge,1)
                myab = find((Gp(BiEdge(k,1),:)+Gp(BiEdge(k,2),:))==2);
                bendlcl(k,:) = [BiEdge(k,:) myab];
            end
            bend(cntb+1:cntb+size(BiEdge,1),:) = Panel{i}(bendlcl);
            cntb = cntb+size(BiEdge,1);
        end
    end
    bend(isnan(bend(:,1)),:)=[];  
    panelctr = [];
else  %'N5B8'
    Nn = size(Node,1);
    bend = nan(6*length(Panel),4);
    panelctr = nan(length(Panel),1);
    cntb = 0;   cntc = 0;
    for i = 1:length(Panel)
        if numel(Panel{i}) == 4
            cntc = cntc+1;
            L1 = norm((Node(Panel{i}(1),:)-Node(Panel{i}(3),:)));
            L2 = norm((Node(Panel{i}(4),:)-Node(Panel{i}(2),:)));
            m = [Node(Panel{i}(3),:)-Node(Panel{i}(1),:)]'/L1;
            n = [Node(Panel{i}(4),:)-Node(Panel{i}(2),:)]'/L2;
            coeff = [m,-n]\[Node(Panel{i}(2),:)-Node(Panel{i}(1),:)]';
            if L1<L2, ctrnode = Node(Panel{i}(1),:)+coeff(1)*m';
            else ctrnode = Node(Panel{i}(2),:)+coeff(2)*n'; end
            % If 1,2,3,4 are co-planar, the two results are identical;
            % if they are not, the central node is placed such that the
            % panel is bended along the shorter diagonal.
            Node(Nn+cntc,:) = ctrnode;
            panelctr(i)=Nn+cntc;
            for k=1:numel(Panel{i})
                cntb = cntb+1;
                ref1 = mod(k-2,numel(Panel{i}))+1;
                ref2 = mod(k,numel(Panel{i}))+1;
                bend(cntb,:) = [Nn+cntc,Panel{i}(k),Panel{i}(ref1),Panel{i}(ref2)];
            end
        elseif numel(Panel{i}) > 4
            cntc = cntc+1;
            ctrnode = FindOptCenter(Node(Panel{i},:));
            Node(Nn+cntc,:) = ctrnode;
            panelctr(i) = Nn+cntc;
            for k=1:numel(Panel{i})
                cntb = cntb+1;
                ref1 = mod(k-2,numel(Panel{i}))+1;
                ref2 = mod(k,numel(Panel{i}))+1;
                bend(cntb,:) = [Nn+cntc,Panel{i}(k),Panel{i}(ref1),Panel{i}(ref2)];
            end
        end
    end
    bend(isnan(bend(:,1)),:)=[];
end
end

%--------------------------------------------------------------------------
function BiEdge = DividePolygon(PolyCoord)
% Generalized N4B5 scheme
if size(PolyCoord,1)<=3
    BiEdge = [];
elseif size(PolyCoord,1)>3
    G = triu(ones(size(PolyCoord,1)),2);
    G(1,end) = 0;
    [I,J] = find(G);
    L2 = sum((PolyCoord(J,:)-PolyCoord(I,:)).^2,2);
    [~,IndMin] = min(L2);
    BiEdge = sort([I(IndMin),J(IndMin)]);
    T1 = [1:BiEdge(1),BiEdge(2):size(PolyCoord,1)];
    T2 = BiEdge(1):BiEdge(2);
    BiEdge = [BiEdge;...
              T1(DividePolygon(PolyCoord(T1,:)));...
              T2(DividePolygon(PolyCoord(T2,:)))];
end
end


function [XC] = FindOptCenter(PolyCoord)
% Generalized N5B8 scheme
G = triu(ones(size(PolyCoord,1)),2); % Generate full graph
G(1,end) = 0;
[I,J] = find(G);
L2 = sum((PolyCoord(J,:)-PolyCoord(I,:)).^2,2);

% Define obj func
objfun = @(XC)sum(sqrt(sum(cross((bsxfun(@minus,PolyCoord(J,:),XC))',(bsxfun(@minus,PolyCoord(I,:),XC))').^2,1))'./(L2.^(1)));

% Set initial guesses
[~,Idmin] = min(L2);
XC01 = (PolyCoord(I(Idmin),:)+PolyCoord(J(Idmin),:))/2;
XC02 = sum(PolyCoord,1)/size(PolyCoord,1);

% Find solution with double starts
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
[XC1,objval1] = fminunc(objfun,XC01,options);
[XC2,objval2] = fminunc(objfun,XC02,options);
XC = XC1*(objval1<=objval2)+XC2*(objval2<objval1);
disp(['MinObj ([MidPoint, CenterPoint]):   ', num2str([objval1, objval2])])
end


function [fold, bdry, Trigl] = findfdbd(Panel,bend)
Nn = max(cellfun(@max,Panel)); 
% triangularization
Panelsize = cellfun(@numel,Panel);
Ptri = cell(sum(Panelsize==3),1);
flg = find(Panelsize==3);
for i = 1:sum(Panelsize==3), Ptri{i} = Panel{flg(i)}; end
Triglraw = [bend(:,[1,2,3]);bend(:,[1,2,4]);cell2mat(Ptri)];
Triglraw = sort(Triglraw,2);
Trigl = unique(Triglraw ,'rows');
% formulate connectivity matrix
Comm = sparse(Nn,size(Trigl,1));
for i=1:size(Trigl,1), Comm(Trigl(i,:),i) = true; end;
% search for fold lines
Ge = Comm'*Comm;
[mf, me] = find(triu(Ge==2)); % triangular meshes that share two common nodes
fold = zeros(length(mf),4);
for i=1:length(mf)
    [link,ia,ib] = intersect(Trigl(mf(i),:),Trigl(me(i),:));
    oftpa = setdiff(1:3,ia);
    oftpb = setdiff(1:3,ib);
    fold(i,:) = [link,Trigl(mf(i),oftpa),Trigl(me(i),oftpb)];
end
fdandbd = sort(fold(:,1:2),2);
onlybd = sort(bend(:,1:2),2);
[~,ibd] = intersect(fdandbd,onlybd,'rows');
fold(ibd,:) = [];

% search for boundaries
Edge = sort([Trigl(:,1) Trigl(:,2); Trigl(:,2) Trigl(:,3); Trigl(:,3) Trigl(:,1)],2);
[u,~,n] = unique(Edge ,'rows');
counts = accumarray(n(:), 1);
bdry = u(counts==1,:);
end


function value = getfieldvalues(options,name,defaultval)
if isempty(options)
     value = defaultval;
     return;
end

try
    value = options.(name);
catch ME
    value = defaultval;
end
end

function [Abarj, kpbj] = GetMaterial(Node,List,indexctr,E,nv,t,Fold,Bend,Bdry)
Fold = sort(Fold(:,1:2),2); Nf = size(Fold,1);
Bend = sort(Bend(:,1:2),2); Nb = size(Bend,1);
Bdry = sort(Bdry,2);        
G = E*t^3/(12*(1-nv^2));

if numel(List)==3
    Pairs = sort([List;List([2:end,1])],1);  Pairs = Pairs';
    Lf = sqrt(sum((Node(Pairs(:,2),:)-Node(Pairs(:,1),:)).^2,2));
    S = 0.5*norm(cross((Node(List(2),:)-Node(List(1),:)),(Node(List(3),:)-Node(List(1),:))));
    A = t*2*S/(sum(Lf)*(1-nv));
    [~,indfd,~] = intersect(Fold,Pairs,'rows');
    if numel(indfd)<3
        [~,inddd,~] = intersect(Bdry,Pairs,'rows');
    else
        inddd = [];
    end
    kpbj = [];
    Abarj = [[indfd+Nb;inddd+Nf+Nb],A*ones(numel(List),1)];
elseif numel(List)==4
    Pairs = sort([List;List([2:end,1])],1);  Pairs = Pairs';
    Lf = sqrt(sum((Node(Pairs(:,2),:)-Node(Pairs(:,1),:)).^2,2));
    [~,indfd,efd] = intersect(Fold,Pairs,'rows');
    Spoke = sort([List' ones(numel(List),1)*indexctr],2); 
    Lb = sqrt(sum(bsxfun(@minus,Node(List,:),Node(indexctr,:)).^2,2));
    [~,indbd,~] = intersect(Bend,Spoke,'rows');
    Dsl = (Lb+Lb([3 4 1 2]));
    kb = (G*(Dsl/t).^(1/3))./Dsl;
    kpbj = [indbd,kb];
    W = sum(Lf([1 3]))/2;  H = sum(Lf([2 4]))/2;
    Af = zeros(4,1);
    Af([1 3]) = t*abs(H^2-nv*W^2)/(2*H*(1-nv^2));
    Af([2 4]) = t*abs(W^2-nv*H^2)/(2*W*(1-nv^2));
    Ab = ones(4,1)*(t*nv*(H^2+W^2)^1.5/(2*H*W*(1-nv^2)));
    if numel(indfd)==numel(List)
        Abarj = [[indbd;indfd+Nb],[Ab;Af]];
    elseif numel(indfd)<numel(List)
        [~,inddd,edd] = intersect(Bdry,Pairs,'rows');
        Abarj = [[indbd;indfd+Nb;inddd+Nf+Nb],[Ab;Af(efd);Af(edd)]];
    end
elseif numel(List)>4
    Pairs = sort([List;List([2:end,1])],1);  Pairs = Pairs';
    Lf = sqrt(sum((Node(Pairs(:,2),:)-Node(Pairs(:,1),:)).^2,2));
    [~,indfd,~] = intersect(Fold,Pairs,'rows');
    Spoke = sort([List' ones(numel(List),1)*indexctr],2); 
    Sa = bsxfun(@minus,Node(List,:),Node(indexctr,:));
    Lb = sqrt(sum(Sa.^2,2));
    Sac = cross(Sa,bsxfun(@minus,Node(List([2:end,1]),:),Node(indexctr,:)));
    S = 0.5*sum(sqrt(sum(Sac.^2,2)));
    A = t*2*S/((sum(Lf)+sum(Lb))*(1-nv));
    [~,indbd,~] = intersect(Bend,Spoke,'rows');
    if numel(indfd)==numel(List)
        Abarj = [[indbd;indfd+Nb],A*ones(2*numel(List),1)];
    elseif numel(indfd)<numel(List)
        [~,inddd,~] = intersect(Bdry,Pairs,'rows');
        Abarj = [[indbd;indfd+Nb;inddd+Nf+Nb],A*ones(2*numel(List),1)];
    end
    Sa = bsxfun(@rdivide,Sa,Lb);
    D = Sa*(Sa');
    [~,Id] = min(D,[],1);
    Dsl = (Lb+Lb(Id));
    kb = (G*(Dsl/t).^(1/3))./Dsl;
    kpbj = [indbd,kb];
end
end

function [B, L] = dirc3d(Node,Ele)
Ne = size(Ele,1); Nn = size(Node,1);
D = [Node(Ele(:,2),1)-Node(Ele(:,1),1), Node(Ele(:,2),2)-Node(Ele(:,1),2),...
    Node(Ele(:,2),3)-Node(Ele(:,1),3)];
L = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);
D = [D(:,1)./L D(:,2)./L D(:,3)./L];
B = sparse(repmat((1:Ne)',1,6),[3*Ele(:,1)-2 3*Ele(:,1)-1 3*Ele(:,1),...
           3*Ele(:,2)-2 3*Ele(:,2)-1 3*Ele(:,2)],[D -D],Ne,3*Nn);
B = -B;
end

