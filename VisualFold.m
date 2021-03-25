function VisualFold(U_his,truss,angles,LF_his,instdof,varargin)

if length(varargin)<1, varargin = {'IntensityMap','off'}; end
options = VFParseOptions(varargin{:});
Node = truss.Node; Trigl = truss.Trigl; 
Panel = angles.Panel;
U_his = [truss.U0,U_his];
if ~isempty(options.CableUSi)
    options.CableUSi = [0*options.CableUSi(:,1),options.CableUSi];
end
if strcmpi(options.IntensityMap,'Vertex')
    if min(options.IntensityData)>=0
        Ng = 110;
        ColCol = jet(Ng);
        VIntensityDataInten = [0*options.IntensityData(:,1),options.IntensityData]/max(max(options.IntensityData));
    else
        cmp = hsv(120);
        cmpr = zeros(10,3); cmpb = cmpr;
        cmpr(:,1) = 0.5:0.05:0.95;
        cmpb(:,3) = 0.95:-0.05:0.5;
        ColCol = [cmpr;cmp(1:21,:);cmp(22:3:61,:);cmp(62:81,:);cmpb];
        Ng = size(ColCol,1);
        VIntensityDataInten = 0.5*(min(max([0*options.IntensityData(:,1),options.IntensityData]/max(max(abs(options.IntensityData))),-1),1)+1);
    end
    VIntensityDataIntenIndex = ceil((Ng-10)*VIntensityDataInten)+9;
elseif strcmpi(options.IntensityMap,'Edge')
    if min(options.IntensityData)>=0
        Ng = 110;
        ColCol = parula(Ng);
        EIntensityDataInten = [0*options.IntensityData(:,1),options.IntensityData]/max(max(options.IntensityData));
    else
        cmp = hsv(120);
        cmpr = zeros(10,3); cmpb = cmpr;
        cmpr(:,1) = 0.5:0.05:0.95;
        cmpb(:,3) = 0.95:-0.05:0.5;
        ColCol = [cmpr;cmp(1:21,:);cmp(22:3:61,:);cmp(62:81,:);cmpb];
        Ng = size(ColCol,1);
        EIntensityDataInten = 0.5*(min(max([0*options.IntensityData(:,1),options.IntensityData]/max(max(abs(options.IntensityData))),-1),1)+1);
    end
        EIntensityDataIntenIndex = ceil((Ng-10)*EIntensityDataInten)+9;
elseif strcmpi(options.IntensityMap,'Face')
    if min(options.IntensityData)>=0
        Ng = 110;
        ColCol = parula(Ng);
        FIntensityDataInten = [0*options.IntensityData(:,1),options.IntensityData]/max(max(options.IntensityData));
    else
        cmp = hsv(120);
        cmpr = zeros(10,3); cmpb = cmpr;
        cmpr(:,1) = 0.5:0.05:0.95;
        cmpb(:,3) = 0.95:-0.05:0.5;
        ColCol  = [cmpr;cmp(1:21,:);cmp(22:3:61,:);cmp(62:81,:);cmpb];
        Ng = size(ColCol,1);
        FIntensityDataInten = 0.5*(min(max([0*options.IntensityData(:,1),options.IntensityData]/max(max(abs(options.IntensityData))),-1),1)+1);
    end
        FIntensityDataIntenIndex = ceil((Ng-10)*FIntensityDataInten)+9;
end

if isempty(LF_his)
    if strcmpi(options.recordtype,'video')
        opengl('hardware')
        writerObj = VideoWriter(options.filename,'MPEG-4');
        writerObj.FrameRate = 24;
        open(writerObj);
    elseif strcmpi(options.recordtype,'imggif')
        filename = [options.filename '.gif'];
    else
        disp('Not recording');
    end
    f1 = figure('units','pixels','position',100+[0 0 720 500]);
    f1.Color = 'w';
    hold on
    for i = 1:size(U_his,2)
        U = U_his(:,i);
        clf
        Nodew = Node;
        Nodew(:,1) = Node(:,1)+U(1:3:end);
        Nodew(:,2) = Node(:,2)+U(2:3:end);
        Nodew(:,3) = Node(:,3)+U(3:3:end);
        
        if strcmpi(options.showinitial,'on')
            PlotOri(Node,Panel,Trigl,'FoldEdgeStyle','-','EdgeShade',0.3,'PanelColor','none');
        end
        if isempty(options.CableUSi)
            if strcmpi(options.IntensityMap,'Vertex')
                PlotOri(Nodew,Panel,Trigl,'FaceVertexColor',ColCol(VIntensityDataIntenIndex(:,i),:))
            elseif strcmpi(options.IntensityMap,'Edge')
                PlotOri(Nodew,Panel,Trigl,'EdgeColor',ColCol(EIntensityDataIntenIndex(:,i),:),'Bars',truss.Bars,...
                        'NumBendHinge',size(angles.bend,1))
            elseif strcmpi(options.IntensityMap,'Face')
                PlotOri(Nodew,Panel,Trigl,'FaceVertexColor',ColCol(FIntensityDataIntenIndex(:,i),:))
            end
        else
            if strcmpi(options.IntensityMap,'Vertex')
                if isfield(truss,'GBars')
                    PlotOri(Nodew,Panel,Trigl,'FaceVertexColor',ColCol(VIntensityDataIntenIndex(:,i),:),'Cables',truss.Cables,'CableSlack',options.CableUSi(:,i),'GBars',truss.GBars)
                else
                    PlotOri(Nodew,Panel,Trigl,'FaceVertexColor',ColCol(VIntensityDataIntenIndex(:,i),:),'Cables',truss.Cables,'CableSlack',options.CableUSi(:,i));
                end
            elseif strcmpi(options.IntensityMap,'Edge')
                PlotOri(Nodew,Panel,Trigl,'EdgeColor',ColCol(EIntensityDataIntenIndex(:,i),:),'Bars',truss.Bars,...
                        'NumBendHinge',size(angles.bend,1),'Cables',truss.Cables,'CableSlack',options.CableUSi(:,i))
            elseif strcmpi(options.IntensityMap,'Face')
                PlotOri(Nodew,Panel,Trigl,'FaceVertexColor',ColCol(FIntensityDataIntenIndex(:,i),:),'Cables',truss.Cables,'CableSlack',options.CableUSi(:,i),'GBars',truss.GBars)
            end
        end

        axis equal; axis off; 
        view(options.viewangle)
        light
        pause(options.pausetime);
        if strcmpi(options.recordtype,'imggif')
            frame = getframe(f1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);                
            if i == 1
                imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', options.pausetime);
            end  
        elseif strcmpi(options.recordtype,'video')
            frame = getframe(f1);
            writeVideo(writerObj,frame);
        end
    end
    hold off;
    if strcmpi(options.recordtype,'video'), close(writerObj); end
    
else
    LF_his = [0*LF_his(1,:);LF_his];
    if size(LF_his,2)>1, LF_his = sum(LF_his,2); end
    if strcmpi(options.recordtype,'video')
        opengl('hardware')
        writerObj = VideoWriter(options.filename,'MPEG-4');
        writerObj.FrameRate = 24;
        open(writerObj);
    elseif strcmpi(options.recordtype,'imggif')
        filename = [options.filename '.gif'];
    else
        disp('Not recording');
    end
    f1 = figure('units','pixels','position',100+[0 0 720 500]);
    f1.Color = 'w';
    hold on
    for i = 1:size(U_his,2)
        U = U_his(:,i);
        clf
        view(options.viewangle)
        Nodew = Node;
        Nodew(:,1) = Node(:,1)+U(1:3:end);
        Nodew(:,2) = Node(:,2)+U(2:3:end);
        Nodew(:,3) = Node(:,3)+U(3:3:end);

        if strcmpi(options.showinitial,'on')
            PlotOri(Node,Panel,Trigl,'FoldEdgeStyle','-','EdgeShade',0.3,'PanelColor','none');
        end
        if strcmpi(options.IntensityMap,'Vertex')
            PlotOri(Nodew,Panel,Trigl,'FaceVertexColor',ColCol(VIntensityDataIntenIndex(:,i),:))
        elseif strcmpi(options.IntensityMap,'Edge')
            PlotOri(Nodew,Panel,Trigl,'EdgeColor',ColCol(EIntensityDataIntenIndex(:,i),:),'Bars',truss.Bars,...
                    'NumBendHinge',size(angles.bend,1))
        elseif strcmpi(options.IntensityMap,'Face')
            PlotOri(Nodew,Panel,Trigl,'FaceVertexColor',ColCol(FIntensityDataIntenIndex(:,i),:))
        else
            PlotOri(Nodew,Panel,Trigl,'Cables',truss.Cables,'CableSlack',options.CableUSi(:,i));
        end
        plot3(Nodew(instdof(1),1),Nodew(instdof(1),2),Nodew(instdof(1),3),'rv','linewidth',2,'MarkerSize',10)
        axis equal; axis off; 
        light
        pause(options.pausetime);
        if strcmpi(options.recordtype,'imggif')
            frame = getframe(f1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);                
            if i == 1
                imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', options.pausetime);
            end  
        elseif strcmpi(options.recordtype,'video')
            frame = getframe(f1);
            writeVideo(writerObj,frame);
        end
    end
    hold off
    if strcmpi(options.recordtype,'video'), close(writerObj); end
    
    if strcmpi(options.recordtype,'video')
        opengl('hardware')
        writerObj = VideoWriter([options.filename '_dispvslambda'],'MPEG-4');
        writerObj.FrameRate = 24;
        open(writerObj);
    elseif strcmpi(options.recordtype,'imggif')
        filename = [options.filename 'dispvslambda' '.gif'];
    end
    f2 = figure('units','pixels','position',100+[0 0 720 500]);
    f2.Color = 'w';
    dsp = sign(instdof(2))*U_his((instdof(1)*3-(3-abs(instdof(2)))),:);
    for i = 1:numel(LF_his)
        clf;
        plot(dsp(1:i),LF_his(1:i),'b-','linewidth',2);
        hold on;
        plot(dsp(i),LF_his(i),'ro','linewidth',2);
        axis(options.axislim);
        xlabel('displacement','fontsize',16)
        ylabel('load factor','fontsize',16)
        pause(options.pausetime)
        if strcmpi(options.recordtype,'imggif')
            frame = getframe(f2);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);                
            if i == 1
                imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', options.pausetime);
            end  
        elseif strcmpi(options.recordtype,'video')
            frame = getframe(f2);
            writeVideo(writerObj,frame);
        end
    end
    hold off
    if strcmp(options.recordtype,'video'), close(writerObj); end
end
end

function options = VFParseOptions(varargin)
IP = inputParser;
IP.addParameter('recordtype', 'none', @ischar)
IP.addParameter('showinitial', 'on', @ischar)
IP.addParameter('filename', 'none', @ischar);
IP.addParameter('pausetime', 0.0001, @isnumeric);
IP.addParameter('axislim',[-inf inf -inf inf], @isnumeric)
IP.addParameter('IntensityMap', 'off', @ischar)
IP.addParameter('IntensityData', [], @isnumeric)
IP.addParameter('viewangle', [35,30], @isnumeric)
IP.addParameter('CableUSi', [], @isnumeric)
IP.parse(varargin{:});
options = IP.Results;
end
    