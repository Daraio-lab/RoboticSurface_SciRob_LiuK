%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 MERLIN2                               %%
%           Written by: Ke Liu (ke.liu@gatech.edu/liuke@caltech.edu)      %
% Ref: K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid    %
%      origami - An efficient computational approach.' PRSA.              %
%      K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to  % 
%      capture highly nonlinear behavior of non-rigid origami.'           %
%      Proceedings of IASS Annual Symposium 2016.                         %
%      E. T. Filipov, K. Liu, T. Tachi, M. Schenk, G. H. Paulino (2017).  %
%      'Bar and hinge models for scalable analysis of origami.'  IJSS     %
%      K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear          %
%      structural analysis of origami assemblages using the MERLIN2       %
%      software.' Origami^7.                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =========== ROBOTIC SURFACE ========================================= %%
clear all; close all; 
%% Define geomtry
% Import geometry in OBJ format
d = 7; Lx = 30; Ly = 30; Nx = 11; Ny = 11; Mx = 4; My = 4; tp = 0.127; ta = 2; wa = 5;% ta = 1.68; wa = 4.27;
[Node, Panel, Cables, GBars] = GetStripe_2D_Mirror(d,Lx,Ly,Nx,Ny,Mx,My,tp,ta,0);
EY = 2.5e3; nu = 0.34; 
Kbb = 5.5*EY*tp^3/(12*(Lx+d)/(Nx));
% Assume that isotropic stretch is 0.7143 (=1.0/1.4)
% groups of LCE strips ("Cables")
XD1 = 1:8;
XD2 = 9:16;
XD3 = 17:24;
YD1 = 25:32;
YD2 = 33:40;
YD3 = 41:48;
XU1 = 49:60;
XU2 = 61:72;
XU3 = 73:84;
YU1 = 85:96;
YU2 = 97:108;
YU3 = 109:120;

% Visualize initial configuration 
figure()
PlotOri(Node+0.1*rand(size(Node)),Panel,[],'Cables',Cables,'Gbars',GBars);
% axis equal; axis off;
light
% Inspect nodal index assignment
figure()
PlotOri(Node,Panel,[],'ShowNumber','on');
% axis equal

%% Set up boundary conditions
Load = [1, 0, 0, 0];
indp = Load(:,1);

% Contact between cable and ribbon
Contact = [];
panelmat = cell2mat(Panel);
for c = 1:size(Cables,1)/2
    ci = Cables(2*c-1,1);
    cj = Cables(2*c-1,2);
    ck = Cables(2*c,1);
    cl = Cables(2*c,2);
    [ip,~] = find(panelmat==ci);
    [jp,~] = find(panelmat==cj);
    [kp,~] = find(panelmat==ck);
    [lp,~] = find(panelmat==cl);
    leftp = intersect(ip,kp);
    rightp = intersect(jp,lp);
    % find ci axis and pointer
    for s = 1:numel(leftp)
        indi = find(panelmat(leftp(s),:)==ci);
        indk = find(panelmat(leftp(s),:)==ck);
        candi = panelmat(leftp(s),mod(indi+(indi-indk)-1,4)+1);
        candk = panelmat(leftp(s),mod(indk+(indk-indi)-1,4)+1);
            
        vec1 = Node(candi,:)-Node(ci,:);
        vec0 = Node(cj,:) - Node(ci,:);
        if (vec1*vec0')>0
            leftcontact = [ci,ck,candi,cj;
                           ck,ci,candk,cl];
            break;
        end
    end
    Contact = [Contact;leftcontact];
    
    for s = 1:numel(rightp)
        indj = find(panelmat(rightp(s),:)==cj);
        indl = find(panelmat(rightp(s),:)==cl);
        candj = panelmat(rightp(s),mod(indj+(indj-indl)-1,4)+1);
        candl = panelmat(rightp(s),mod(indl+(indl-indj)-1,4)+1);
            
        vec1 = Node(candj,:)-Node(cj,:);
        vec0 = Node(ci,:) - Node(cj,:);
        if (vec1*vec0')>0
            rightcontact = [cj,cl,candj,ci;
                            cl,cj,candl,ck];
            break;
        end
    end
    Contact = [Contact;rightcontact];
end

Kc = 1e-5+zeros(size(Contact,1),1);
pc0 = zeros(size(Contact,1),1); 
for i = 1:size(pc0,1)
    pci = FoldKe(Node,Contact(i,:));
    if pci>=pi
        pc0(i) = 2*pi-pi/6; 
    else
        pc0(i) = pi/6; 
    end
end

ShrinkAssign = zeros(size(Cables,1),1);
% Saddle
Ind_1 = 132;
Ind_2 = 145;
Ind_3 = 144;
Ind_4 = 157;
Supp = [ Ind_1, 1, 1, 1;
         Ind_2, 1, 1, 1;
         Ind_3, 1, 1, 1;
         Ind_4, 1, 1, 1];
ShrinkAssign(XD1) = 0.75;
ShrinkAssign(XD2) = 0.75;
ShrinkAssign(XD3) = 0.75;
ShrinkAssign(YD1) = 0.95;
ShrinkAssign(YD2) = 0.95;
ShrinkAssign(YD3) = 0.95;
ShrinkAssign(XU1) = 0.95;
ShrinkAssign(XU2) = 0.95;
ShrinkAssign(XU3) = 0.95;
ShrinkAssign(YU1) = 0.75;
ShrinkAssign(YU2) = 0.75;
ShrinkAssign(YU3) = 0.75;

% Dome 
% [~,Ind_1] = min(sum((Node-[30,0,0]).^2,2));
% [~,Ind_2] = min(sum((Node-[30,141,0]).^2,2));
% [~,Ind_3] = min(sum((Node-[111,0,0]).^2,2));
% [~,Ind_4] = min(sum((Node-[111,141,0]).^2,2));
% [~,Ind_5] = min(sum((Node-[141,30,0]).^2,2));
% [~,Ind_6] = min(sum((Node-[141,111,0]).^2,2));
% [~,Ind_7] = min(sum((Node-[0,30,0]).^2,2));
% [~,Ind_8] = min(sum((Node-[0,111,0]).^2,2));
% Supp = [ Ind_1, 1, 1, 1;
%          Ind_2, 1, 0, 1;
%          Ind_3, 0, 1, 1;
%          Ind_4, 0, 0, 1;
%          Ind_5, 0, 0, 1;
%          Ind_6, 0, 0, 1;
%          Ind_7, 0, 0, 1;
%          Ind_8, 0, 0, 1];
% ShrinkAssign(XD1) = 0.75;
% ShrinkAssign(XD2) = 0.75;
% ShrinkAssign(XD3) = 0.75;
% ShrinkAssign(YD1) = 0.75;
% ShrinkAssign(YD2) = 0.75;
% ShrinkAssign(YD3) = 0.75;
% ShrinkAssign(XU1) = 0.95;
% ShrinkAssign(XU2) = 0.95;
% ShrinkAssign(XU3) = 0.95;
% ShrinkAssign(YU1) = 0.95;
% ShrinkAssign(YU2) = 0.95;
% ShrinkAssign(YU3) = 0.95;

% Cynlinder
% [~,Ind_1] = min(sum((Node-[30,0,0]).^2,2));
% [~,Ind_2] = min(sum((Node-[30,141,0]).^2,2));
% [~,Ind_3] = min(sum((Node-[111,0,0]).^2,2));
% [~,Ind_4] = min(sum((Node-[111,141,0]).^2,2));
% [~,Ind_5] = min(sum((Node-[141,30,0]).^2,2));
% [~,Ind_6] = min(sum((Node-[141,111,0]).^2,2));
% [~,Ind_7] = min(sum((Node-[0,30,0]).^2,2));
% [~,Ind_8] = min(sum((Node-[0,111,0]).^2,2));
% Supp = [ Ind_1, 1, 1, 1;
%          Ind_2, 1, 0, 1;
%          Ind_3, 0, 1, 1;
%          Ind_4, 0, 0, 1;
%          Ind_5, 0, 0, 1;
%          Ind_6, 0, 0, 1;
%          Ind_7, 0, 0, 1;
%          Ind_8, 0, 0, 1];
% ShrinkAssign(XD1) = 0.95;
% ShrinkAssign(XD2) = 0.95;
% ShrinkAssign(XD3) = 0.95;
% ShrinkAssign(YD1) = 0.95;
% ShrinkAssign(YD2) = 0.95;
% ShrinkAssign(YD3) = 0.95;
% ShrinkAssign(XU1) = 0.95;
% ShrinkAssign(XU2) = 0.95;
% ShrinkAssign(XU3) = 0.95;
% ShrinkAssign(YU1) = 0.75;
% ShrinkAssign(YU2) = 0.75;
% ShrinkAssign(YU3) = 0.75;

% Table
% [~,Ind_1] = min(sum((Node-[30,0,0]).^2,2));
% [~,Ind_2] = min(sum((Node-[30,141,0]).^2,2));
% [~,Ind_3] = min(sum((Node-[111,0,0]).^2,2));
% [~,Ind_4] = min(sum((Node-[111,141,0]).^2,2));
% [~,Ind_5] = min(sum((Node-[141,30,0]).^2,2));
% [~,Ind_6] = min(sum((Node-[141,111,0]).^2,2));
% [~,Ind_7] = min(sum((Node-[0,30,0]).^2,2));
% [~,Ind_8] = min(sum((Node-[0,111,0]).^2,2));
% Supp = [ Ind_1, 1, 1, 1;
%          Ind_2, 1, 0, 1;
%          Ind_3, 0, 1, 1;
%          Ind_4, 0, 0, 1;
%          Ind_5, 0, 0, 1;
%          Ind_6, 0, 0, 1;
%          Ind_7, 0, 0, 1;
%          Ind_8, 0, 0, 1];
% ShrinkAssign(XD1) = 0.75;
% ShrinkAssign(XD2) = 0.75;
% ShrinkAssign(XD3) = 0.75;
% ShrinkAssign(YD1) = 0.75;
% ShrinkAssign(YD2) = 0.75;
% ShrinkAssign(YD3) = 0.75;
% ShrinkAssign(XU1) = 0.72; % to compensate the GBars
% ShrinkAssign(XU2) = 0.72;
% ShrinkAssign(XU3) = 0.72;
% ShrinkAssign(YU1) = 0.72;
% ShrinkAssign(YU2) = 0.72;
% ShrinkAssign(YU3) = 0.72;

% Cone
% Ind_1 = 132;
% Ind_2 = 145;
% Ind_3 = 144;
% Ind_4 = 157;
% Supp = [ Ind_1, 1, 1, 1;
%          Ind_2, 1, 1, 1;
%          Ind_3, 1, 1, 1;
%          Ind_4, 1, 1, 1];
% ShrinkAssign(XD1) = 0.75;
% ShrinkAssign(XD2) = 0.85;
% ShrinkAssign(XD3) = 0.95;
% ShrinkAssign(YD1) = 0.9;
% ShrinkAssign(YD2) = 0.9;
% ShrinkAssign(YD3) = 0.9;
% ShrinkAssign(XU1) = 0.95;
% ShrinkAssign(XU2) = 0.95;
% ShrinkAssign(XU3) = 0.95;
% ShrinkAssign(YU1) = 0.9;
% ShrinkAssign(YU2) = 0.9;
% ShrinkAssign(YU3) = 0.9;

% Dustpan
% Ind_1 = 132;
% Ind_2 = 145;
% Ind_3 = 144;
% Ind_4 = 157;
% Supp = [ Ind_1, 1, 1, 1;
%          Ind_2, 1, 1, 1;
%          Ind_3, 1, 1, 1;
%          Ind_4, 1, 1, 1];
% ShrinkAssign(XD1([1:4])) = 0.95;
% ShrinkAssign(XD1([5:8])) = 0.85;
% ShrinkAssign(XD2([1:4])) = 0.95;
% ShrinkAssign(XD2([5:8])) = 0.85;
% ShrinkAssign(XD3([1:4])) = 0.95;
% ShrinkAssign(XD3([5:8])) = 0.85; 
% ShrinkAssign(XU1([1:4,5,7])) = 0.8;
% ShrinkAssign(XU1([6,8,9:12])) = 0.95;
% ShrinkAssign(XU2([1:4,5,7])) = 0.8;
% ShrinkAssign(XU2([6,8,9:12])) = 0.95;
% ShrinkAssign(XU3([1:4,5,7])) = 0.8;
% ShrinkAssign(XU3([6,8,9:12])) = 0.95;
% ShrinkAssign(YD1) = 0.95;
% ShrinkAssign(YD2) = 0.95;
% ShrinkAssign(YD3) = 0.95;
% ShrinkAssign(YU1) = 0.8;
% ShrinkAssign(YU2) = 0.8;
% ShrinkAssign(YU3) = 0.8;

% Wave
% Ind_1 = 132;
% Ind_2 = 145;
% Ind_3 = 144;
% Ind_4 = 157;
% Supp = [ Ind_1, 1, 1, 1;
%          Ind_2, 1, 1, 1;
%          Ind_3, 1, 1, 1;
%          Ind_4, 1, 1, 1];
% ShrinkAssign(XD1([1:4])) = 0.95;
% ShrinkAssign(XD1([5:8])) = 0.75;
% ShrinkAssign(XD2([1:4])) = 0.95;
% ShrinkAssign(XD2([5:8])) = 0.75;
% ShrinkAssign(XD3([1:4])) = 0.95;
% ShrinkAssign(XD3([5:8])) = 0.75;
% ShrinkAssign(XU1([1:4,5,7])) = 0.75;
% ShrinkAssign(XU1([6,8,9:12])) = 0.95;
% ShrinkAssign(XU2([1:4,5,7])) = 0.75;
% ShrinkAssign(XU2([6,8,9:12])) = 0.95;
% ShrinkAssign(XU3([1:4,5,7])) = 0.75;
% ShrinkAssign(XU3([6,8,9:12])) = 0.95;
% ShrinkAssign(YD1) = 0.9;
% ShrinkAssign(YD2) = 0.9;
% ShrinkAssign(YD3) = 0.9;
% ShrinkAssign(YU1) = 0.9;
% ShrinkAssign(YU2) = 0.9;
% ShrinkAssign(YU3) = 0.9;

% Twisted
% Ind_1 = 132;
% Ind_2 = 145;
% Ind_3 = 144;
% Ind_4 = 157;
% Supp = [ Ind_1, 1, 1, 1;
%          Ind_2, 1, 1, 1;
%          Ind_3, 1, 1, 1;
%          Ind_4, 1, 1, 1];
% ShrinkAssign(XD1) = 0.85;
% ShrinkAssign(XD2) = 0.8;
% ShrinkAssign(XD3) = 0.95;
% ShrinkAssign(YD1) = 0.9;
% ShrinkAssign(YD2) = 0.9;
% ShrinkAssign(YD3) = 0.9;
% ShrinkAssign(XU1) = 0.95;
% ShrinkAssign(XU2) = 0.8;
% ShrinkAssign(XU3) = 0.85;
% ShrinkAssign(YU1) = 0.9;
% ShrinkAssign(YU2) = 0.9;
% ShrinkAssign(YU3) = 0.9;


%% Define material and modeling parameters
% Simulation options using the N5B8 model
AnalyInputOpt = struct(...
    'ModelType','N5B8',...
    'MaterCalib','auto',...  
    'ModElastic', EY,...
    'Poisson', nu,...
    'Thickness', tp,... 
    'Kf', Kbb, ...
    'Kb', Kbb, ...
    'RotSprBend', @(he,h0,Kb,L0)EnhancedLinear(he,h0,Kb,L0,30,330),...
    'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,30,330),...
    'LoadType','Displacement',...    % Displacement load
    'DispStep',300, ... 
    'Cables', Cables,...
    'CableCM', @(Ex, Pstr)OgdenLCE(Ex, Pstr), ...
    'CableA', wa*ta/2,...
    'GBars', GBars,...
    'GBarCM', @(Ex, Pstr)Ogden(Ex, 5e3*2),...
    'GBarA', d*tp,...
    'TargetPrestrain', 1-ShrinkAssign,...
    'MaxIcr', 2000);

%% Perform analysis
% Assemble input data
[truss, angles, AnalyInputOpt] = PrepareData(Node,Panel,Supp,Load,AnalyInputOpt);
angles.bend = [angles.bend;Contact];
angles.Kb = [angles.Kb;Kc];
angles.pb0 = [angles.pb0;pc0];
% Specify initial deformation state
truss.U0 = zeros(3*size(truss.Node,1),1);
% Perform path-following analysis
[Uhis,Fhis,truss] = PathAnalysis(truss,angles,AnalyInputOpt);
% Postprocess output data
Uhis = real(Uhis);
Fhis = real(Fhis);
STAT = PostProcess(Uhis,truss,angles); 

%% Visualize simulation
instdof = [indp(1),1];
interv = 1; endicrm = size(Uhis,2);
% Animation monitoring node-wise change
VIntensityDataInten = zeros(size(truss.Node,1),size(Uhis,2));
IntensityDataM = bsxfun(@times,STAT.bar.Sx,truss.A);
for k = 1:size(Uhis,2)
    IntensityDataIntenk = sparse(truss.Bars(:,1),truss.Bars(:,2),abs(IntensityDataM(:,k)),size(truss.Node,1),size(truss.Node,1));
    VIntensityDataInten(:,k) = sum((IntensityDataIntenk+IntensityDataIntenk'),2); 
end
VisualFold(Uhis(:,1:interv:endicrm),truss,angles,[],[],'IntensityMap','Vertex','IntensityData',VIntensityDataInten,'CableUSi',-STAT.cable.Ex)
