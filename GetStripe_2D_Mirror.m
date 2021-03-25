function [NODE, PANEL, CABLE, GBAR] = GetStripe_2D_Mirror(d,Lx,Ly,Nx,Ny,Mx,My,tp,ta,AlternateFlag)

gap = tp+ta/2;

%Nx, Ny = Number of Panels within a segments
%Mx, My = Number of segments in x, y direction 
%Lx, Ly = Length between 2 stripes in x, y direction

% total number nodes, panels, cables
%nodes_tot = 2*stripes_x*(px+1)+2*(My-1)*stripes_y*Ny; 
%panels_tot = px*stripes_x + My*stripes_y;
%cable_tot = (stripes_x-1)*stripes_y+stripes_x*(stripes_y-1);

NODE = zeros(2*Mx*(Nx+1),3);
TotalRibbonNodes = 2*Mx*(Nx+1)*(My-1)+2*My*(Ny+1)*(Mx-1)-4*(My-1)*(Mx-1);
NODEJ = zeros((Mx-1)*(My-1)*4,3);
PANEL = cell(Mx*Nx+Mx-1,1);
CABLE = zeros(2*Mx,2);
CABLE_UPX = zeros(2*2*(Mx-1),2);
CABLE_UPX_TEMP = zeros(2*2*(Mx-1),2);
CABLE_UPY = zeros(2*2*(My-1)*(Mx-1),2);
GBAR = zeros(6*(Mx-1),2);

% CREATE FIRST STRIPE IN X-DIRECTION

if mod(Nx+1,2)~=0, error('Upper strings cannot be attached!'); end

spacings = [linspace(0,(Lx-d)/2,(Nx+1)/2),linspace((Lx+d)/2,Lx,(Nx+1)/2)];

for i=1:Mx
    culm = (i-1)*2*(Nx+1);
    NODE(culm+(1:Nx+1),1) = spacings+(i-1)*(Lx+d);
    NODE(culm+((Nx+2):(2*(Nx+1))),1) = spacings+(i-1)*(Lx+d); 
    NODE(culm+((Nx+2):(2*(Nx+1))),2) = d;
    NODE(culm+(1:Nx+1),3) = gap*sin(pi*spacings/Lx).^2;
    NODE(culm+((Nx+2):(2*(Nx+1))),3) = gap*sin(pi*spacings/Lx).^2;

    for j=1:Nx
        PANEL{(i-1)*Nx+j} = [j j+Nx+1 j+Nx+2 j+1]+culm;
    end
    if i>1
        PANEL{Nx*Mx+i-1} = [culm-Nx-1 culm culm+Nx+2 culm+1];
    end
    
    CABLE((i-1)*2+1,:) = [culm+1, culm+Nx+1];
    CABLE((i-1)*2+2,:) = [culm+1, culm+Nx+1]+Nx+1;
    if i>1
        NODEJ((i-2)*4+(1:4),:) = NODE([culm-Nx-1 culm culm+Nx+2 culm+1],:);
        NODEJ((i-2)*4+(1:4),3) = NODE((Nx+1)/2+1+culm-(Nx+1),3);
        CABLE_UPX((i-2)*4+(1:2),:) = [(Nx+1)/2+1+culm-2*(Nx+1), TotalRibbonNodes+(i-2)*4+1;           
                                    TotalRibbonNodes+(i-2)*4+4, (Nx+1)/2+culm];
        CABLE_UPX((i-2)*4+(3:4),:) = [  (Nx+1)/2+1+culm-(Nx+1), TotalRibbonNodes+(i-2)*4+2;        
                                    TotalRibbonNodes+(i-2)*4+3, (Nx+1)/2+culm+Nx+1];
        GBAR((i-2)*6+(1:6),:) = TotalRibbonNodes+(i-2)*4+[1 3; 2 4; 1 2; 2 3; 3 4; 4 1]; 
    end
end

% SHIFT X STRIPES IN Y DIRECTION

NODE(:,2) = NODE(:,2)+Ly;
NODEJ(1:4*(Mx-1),2) = NODEJ(1:4*(Mx-1),2)+Ly;
NODE_TEMP = NODE;
PANEL_TEMP = PANEL;

for l=1:My-2
    
   NODE_SHIFTED = NODE_TEMP;
   NODE_SHIFTED(:,2) = NODE_TEMP(:,2)+l*(Ly+d);
   NODEJ(4*l*(Mx-1)+(1:4*(Mx-1)),:) = NODEJ(1:4*(Mx-1),:);
   NODEJ(4*l*(Mx-1)+(1:4*(Mx-1)),2) = NODEJ(1:4*(Mx-1),2)+l*(Ly+d);
  
   NODE = [NODE; NODE_SHIFTED];
   
   PANEL_SHIFTED = PANEL_TEMP;

    for p=1:(Mx*Nx+Mx-1)
        PANEL_SHIFTED{p}=PANEL_SHIFTED{p}+l*size(NODE_TEMP,1);
    end

    PANEL = [PANEL;PANEL_SHIFTED];
    
    CABLE = [CABLE;CABLE(1:2*Mx,:)+l*size(NODE_TEMP,1)];
    
    for lx = 1:Mx-1
        CABLE_UPX_TEMP((lx-1)*4+(1:4),:) = CABLE_UPX((lx-1)*4+(1:4),:)+l*(Mx-1)*4;
        CABLE_UPX_TEMP((lx-1)*4+1,1) = CABLE_UPX((lx-1)*4+1,1)+l*size(NODE_TEMP,1);
        CABLE_UPX_TEMP((lx-1)*4+2,2) = CABLE_UPX((lx-1)*4+2,2)+l*size(NODE_TEMP,1);
        CABLE_UPX_TEMP((lx-1)*4+3,1) = CABLE_UPX((lx-1)*4+3,1)+l*size(NODE_TEMP,1);
        CABLE_UPX_TEMP((lx-1)*4+4,2) = CABLE_UPX((lx-1)*4+4,2)+l*size(NODE_TEMP,1);    
    end
    
    CABLE_UPX = [CABLE_UPX;CABLE_UPX_TEMP];
    
    GBAR = [GBAR;GBAR(1:6*(Mx-1),:)+l*4*(Mx-1)];
end

% k = 2*(Mx*(Nx+1)+2)*(My+1);
k = size(NODE,1);

%% Y
spacings_y =  [linspace(0,(Ly-d)/2,(Ny+1)/2),linspace((Ly+d)/2,Ly,(Ny+1)/2)];
for ix=1:Mx-1
    shift_x = ix*(Lx+d)-d;
    for iy=1:My
        if iy == 1
            NODE_Y_Bottom = [shift_x,0,0;shift_x+d,0,0];
            NODE = [NODE;NODE_Y_Bottom];
        elseif iy == My
            NODE_Y_Top = [shift_x,(My*Ly+(My-1)*d),0;shift_x+d,(My*Ly+(My-1)*d),0];
            NODE = [NODE;NODE_Y_Top];
        end
        NODE_Y_SEG = zeros(2*(Ny-1),3);
        shift_y = (iy-1)*(Ly+d);
        NODE_Y_SEG(:,2) = [spacings_y(2:end-1),spacings_y(2:end-1)];
        NODE_Y_SEG(Ny:2*(Ny-1),1) = d;
        if AlternateFlag
            NODE_Y_SEG(1:Ny-1,3) = -gap*sin(pi*spacings_y(2:end-1)/Ly).^2;
            NODE_Y_SEG(Ny:2*(Ny-1),3) = -gap*sin(pi*spacings_y(2:end-1)/Ly).^2;
        else
            NODE_Y_SEG(1:Ny-1,3) = gap*sin(pi*spacings_y(2:end-1)/Ly).^2;
            NODE_Y_SEG(Ny:2*(Ny-1),3) = gap*sin(pi*spacings_y(2:end-1)/Ly).^2;
        end
        NODE_Y_SEG(:,1) = NODE_Y_SEG(:,1)+shift_x;
        NODE_Y_SEG(:,2) = NODE_Y_SEG(:,2)+shift_y;

        [~,IND_NODE_ll] = min(abs(NODE(:,1)-NODE_Y_SEG(1,1))+abs(NODE(:,2)-NODE_Y_SEG(1,2)));
        [~,IND_NODE_lr] = min(abs(NODE(:,1)-NODE_Y_SEG(Ny,1))+abs(NODE(:,2)-NODE_Y_SEG(Ny,2)));
        [~,IND_NODE_ul] = min(abs(NODE(:,1)-NODE_Y_SEG(Ny-1,1))+abs(NODE(:,2)-NODE_Y_SEG(Ny-1,2)));
        [~,IND_NODE_ur] = min(abs(NODE(:,1)-NODE_Y_SEG(2*(Ny-1),1))+abs(NODE(:,2)-NODE_Y_SEG(2*(Ny-1),2)));
        PANEL_Y_SEG = cell(Ny,1);
        PANEL_Y_SEG{1} = [IND_NODE_ll,IND_NODE_lr,size(NODE,1)+Ny,size(NODE,1)+1];
        PANEL_Y_SEG{Ny} = [size(NODE,1)+Ny-1,size(NODE,1)+2*(Ny-1),IND_NODE_ur,IND_NODE_ul];
        for j=1:Ny-2
            PANEL_Y_SEG{j+1} = [j, j+(Ny-1), j+(Ny-1)+1, j+1]+size(NODE,1);
        end
        PANEL = [PANEL; PANEL_Y_SEG];
        
        CABLE_Y_SEG_L = [IND_NODE_ll,IND_NODE_ul;
                         IND_NODE_lr,IND_NODE_ur];
        CABLE = [CABLE;CABLE_Y_SEG_L];
        if iy>1 
            if iy == 2
                CABLE_UPY(4*(My-1)*(ix-1)+4*(iy-2)+(1:2),:) = [(Ny-1)/2+size(NODE,1), TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+2;
%                                          TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+2, TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+1;   
                                         TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+1, (Ny-1)/2-2*(Ny-1)+1+size(NODE,1)];
                CABLE_UPY(4*(My-1)*(ix-1)+4*(iy-2)+(3:4),:) = [(Ny-1)/2+Ny-1+size(NODE,1), TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+3;
%                                               TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+3, TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+4;
                                              TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+4, (Ny-1)/2-(Ny-1)+1+size(NODE,1)];
            elseif iy == My
                CABLE_UPY(4*(My-1)*(ix-1)+4*(iy-2)+(1:2),:) = [(Ny-1)/2+size(NODE,1), TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+2;
%                                          TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+2, TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+1;                                            
                                         TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+1, (Ny-1)/2-2*(Ny-1)-1+size(NODE,1)];
                CABLE_UPY(4*(My-1)*(ix-1)+4*(iy-2)+(3:4),:) = [(Ny-1)/2+Ny-1+size(NODE,1), TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+3;
%                                               TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+3, TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+4;                                  
                                              TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+4, (Ny-1)/2-(Ny-1)-1+size(NODE,1)];
            else 
                CABLE_UPY(4*(My-1)*(ix-1)+4*(iy-2)+(1:2),:) = [(Ny-1)/2+size(NODE,1), TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+2;
%                                          TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+2, TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+1;   
                                         TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+1, (Ny-1)/2+1-2*(Ny-1)+size(NODE,1)];
                CABLE_UPY(4*(My-1)*(ix-1)+4*(iy-2)+(3:4),:) = [(Ny-1)/2+Ny-1+size(NODE,1), TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+3;
%                                               TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+3, TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+4;                                  
                                              TotalRibbonNodes+(iy-2)*(Mx-1)*4+(ix-1)*4+4, (Ny-1)/2+1-(Ny-1)+size(NODE,1)];
            end
        end
        NODE = [NODE; NODE_Y_SEG];
    end
end

NODE = [NODE;NODEJ];
CABLE = [CABLE; CABLE_UPX; CABLE_UPY];





