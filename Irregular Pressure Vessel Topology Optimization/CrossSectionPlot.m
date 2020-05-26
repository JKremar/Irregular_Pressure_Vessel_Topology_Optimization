clear all
close all
clc
addpath([pwd,'\IrregularShapeSubfunctions'])
%Plots Cross-Sections of Structures

     
folder='0.45finish23-May-2020 10,08,20';    %iIPV 30 v=0.45
iteration='end';    %iteration # or 'end' for last iteration



global Domain W Esize Done
if(ischar(iteration))
    mat_files=dir([folder,'/*.mat']);
    %iteration=numel(mat_files)-1;
    iteration=142;
end
Done=0;
load([folder,'/Iteration0'])
Esize=max(nodes(elements(1,:),:))-min(nodes(elements(1,:),:));
volTot=prod(Esize)*size(elements,1);
load([folder,'/Iteration',num2str(iteration)])
Domain=[min(nodes);max(nodes)];
W=Domain;
Esize=max(nodes(elements(1,:),:))-min(nodes(elements(1,:),:));
[Ke,B,C]=stiff3D(29.5*10^6,0.29,Esize);
Yield=150*10^3; %psi
VonMises=zeros(size(elements,1),1);
StrucSize=size(struc);

ind=find(struc(map));
for(e=1:sum(struc(map)))
    stress=C*B*U(3*repelem(elements(ind(e),:),1,3)-repmat([2,1,0],1,8));
    VonMises(ind(e))=sqrt(sum((stress(1:3)-stress([2,3,1])).^2)+6*sum(stress(4:6).^2))/sqrt(2);
end
CalcComp=~exist('compE','var');
OldNodes=nodes;
Cent=nodes(elements(:,1),:)+Esize/2;
use=find(Cent(:,1)>=W(1,1)&Cent(:,1)<=W(2,1)&...
    Cent(:,2)>=W(1,2)&Cent(:,2)<=W(2,2)&...
    Cent(:,3)>=W(1,3)&Cent(:,3)<=W(2,3));
f=figure('Units','normalized','color','w');
fig=plotstructure(elements(use,:),nodes,struc,map(use));
%fig=plottrans(elements,nodes,struc,map,boundary);
axis equal; axis tight;  view([30,30]);   drawnow;
xlabel('x');    ylabel('y');    zlabel('z');
lgd=legend('Solid');
lgd.Position=[0.85,0.85,0.1,0.1];

%X Limits Control Panel----------------------------------------------------
Xrange=uipanel('Title','X Limits','Position',[0.01,0.775,0.18755,0.125]);
uicontrol(Xrange,'Style','text','String','Max:','Units','normalized','FontUnits','normalized',...
    'Position',[0,0.6,0.185,0.225],'FontSize',0.9);
Xmax_down=uicontrol(Xrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.2,0.6,0.25,0.225],'String','Down','FontSize',0.9,'Callback',@XmaxDPushed);
Xmax_up=uicontrol(Xrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.5,0.6,0.25,0.225],'String','Up','FontSize',0.9,'Callback',@XmaxUPushed);
Xmax_Text=uicontrol(Xrange,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.775,0.6,0.22,0.225],'FontSize',0.9);
uicontrol(Xrange,'Style','text','String','Min:','Units','normalized','FontUnits','normalized',...
    'Position',[0,0.3,0.185,0.225],'FontSize',0.9);
Xmin_down=uicontrol(Xrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.2,0.3,0.25,0.225],'String','Down','FontSize',0.9,'Callback',@XminDPushed);
Xmin_up=uicontrol(Xrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.5,0.3,0.25,0.225],'String','Up','FontSize',0.9,'Callback',@XminUPushed);
Xmin_Text=uicontrol(Xrange,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.775,0.3,0.22,0.225],'FontSize',0.9);
%--------------------------------------------------------------------------

%Y Limits Control Panel----------------------------------------------------
Yrange=uipanel('Title','Y Limits','Position',[0.01,0.625,0.18755,0.125]);
uicontrol(Yrange,'Style','text','String','Max:','Units','normalized','FontUnits','normalized',...
    'Position',[0,0.6,0.185,0.225],'FontSize',0.9);
Ymax_down=uicontrol(Yrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.2,0.6,0.25,0.225],'String','Down','FontSize',0.9,'Callback',@YmaxDPushed);
Ymax_up=uicontrol(Yrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.5,0.6,0.25,0.225],'String','Up','FontSize',0.9,'Callback',@YmaxUPushed);
Ymax_Text=uicontrol(Yrange,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.775,0.6,0.22,0.225],'FontSize',0.9);
uicontrol(Yrange,'Style','text','String','Min:','Units','normalized','FontUnits','normalized',...
    'Position',[0,0.3,0.185,0.225],'FontSize',0.9);
Ymin_down=uicontrol(Yrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.2,0.3,0.25,0.225],'String','Down','FontSize',0.9,'Callback',@YminDPushed);
Ymin_up=uicontrol(Yrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.5,0.3,0.25,0.225],'String','Up','FontSize',0.9,'Callback',@YminUPushed);
Ymin_Text=uicontrol(Yrange,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.775,0.3,0.22,0.225],'FontSize',0.9);
%--------------------------------------------------------------------------

%Z Limits Control Panel----------------------------------------------------
Zrange=uipanel('Title','Z Limits','Position',[0.01,0.475,0.18755,0.125]);
uicontrol(Zrange,'Style','text','String','Max:','Units','normalized','FontUnits','normalized',...
    'Position',[0,0.6,0.185,0.225],'FontSize',0.9);
Zmax_down=uicontrol(Zrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.2,0.6,0.25,0.225],'String','Down','FontSize',0.9,'Callback',@ZmaxDPushed);
Zmax_up=uicontrol(Zrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.5,0.6,0.25,0.225],'String','Up','FontSize',0.9,'Callback',@ZmaxUPushed);
Zmax_Text=uicontrol(Zrange,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.775,0.6,0.22,0.225],'FontSize',0.9);
uicontrol(Zrange,'Style','text','String','Min:','Units','normalized','FontUnits','normalized',...
    'Position',[0,0.3,0.185,0.225],'FontSize',0.9);
Zmin_down=uicontrol(Zrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.2,0.3,0.25,0.225],'String','Down','FontSize',0.9,'Callback',@ZminDPushed);
Zmin_up=uicontrol(Zrange,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.5,0.3,0.25,0.225],'String','Up','FontSize',0.9,'Callback',@ZminUPushed);
Zmin_Text=uicontrol(Zrange,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.775,0.3,0.22,0.225],'FontSize',0.9);
%--------------------------------------------------------------------------

%Deflection Controls-------------------------------------------------------
DeflecTB=uicontrol(f,'Style','togglebutton','Units','normalized','Position',[0.01,0.42,0.1875,0.04],'String','Deflection','FontUnits','normalized','FontSize',0.75);
MagSlide=uicontrol(f,'Style','slider','Units','normalized','Position',[0.01,0.375,0.11,0.04],'Min',0,'Max',100,'Value',10,'SliderStep',[1/1000,0.01]);
MagText=uicontrol(f,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.12,0.375,0.08,0.04],'FontSize',0.75);
set(MagText,'String',sprintf('Mag:%3.1f',MagSlide.Value))
%--------------------------------------------------------------------------
%Plot Type Controls--------------------------------------------------------
P=uibuttongroup(f,'Position',[0.01,0.135,0.18755,0.21875],'Units','normalized');%'SelectionChangedFcn',@Ptype
tb1=uicontrol(P,'Style','togglebutton','Units','normalized','Position',[0.05,0.76,0.9,0.2],'String','Solid','FontUnits','normalized','FontSize',0.75);
tb2=uicontrol(P,'Style','togglebutton','Units','normalized','Position',[0.05,0.52,0.9,0.2],'String','Transparent','FontUnits','normalized','FontSize',0.75);
tb3=uicontrol(P,'Style','togglebutton','Units','normalized','Position',[0.05,0.28,0.9,0.2],'String','Stress','FontUnits','normalized','FontSize',0.75);
tb4=uicontrol(P,'Style','togglebutton','Units','normalized','Position',[0.05,0.04,0.9,0.2],'String','Void Only','FontUnits','normalized','FontSize',0.75);
%--------------------------------------------------------------------------

S=uicontrol(f,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.01,0.0775,0.1875,0.05],'String','Save','FontSize',0.9,'Callback',@SPushed);
D=uicontrol(f,'Style','pushbutton','Units','normalized','FontUnits','normalized',...
    'Position',[0.01,0.015,0.1875,0.05],'String','Done','FontSize',0.9,'Callback',@DPushed);
uicontrol(f,'Style','text','String','Iteration:','Units','normalized','FontUnits','normalized',...
    'Position',[0.025,0.94,0.15,0.05],'FontSize',0.9);
itr=uicontrol(f,'Style','edit','Units','normalized','FontUnits','normalized',...
    'Position',[0.18,0.94,0.0575,0.05],'String',num2str(iteration),'FontSize',0.9);

%Objective and Constraint Values-------------------------------------------
CompTot=uicontrol(f,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.7375,0.15,0.25,0.025],'FontSize',0.9,'horizontalAlignment','right');
CompWin=uicontrol(f,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.7375,0.12,0.25,0.025],'FontSize',0.9,'horizontalAlignment','right');
VolTot=uicontrol(f,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.7375,0.09,0.25,0.025],'FontSize',0.9,'horizontalAlignment','right');
Volfrac=uicontrol(f,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.7375,0.06,0.25,0.025],'FontSize',0.9,'horizontalAlignment','right');
Intfrac=uicontrol(f,'Style','text','Units','normalized','FontUnits','normalized',...
    'Position',[0.7375,0.03,0.25,0.025],'FontSize',0.9,'horizontalAlignment','right');
set(VolTot,'String',sprintf('Total Volume Fraction:%1.3f',prod(Esize)*sum(struc(map))/volTot))
set(Volfrac,'String',sprintf('Window Volume Fraction:%1.3f',prod(Esize)*sum(struc(map(use)))/volTot))
set(Intfrac,'String',sprintf('Window Interior Volume Fraction:%1.3f',prod(Esize)*sum(struc(map(setdiff(use,boundary))))/volTot))
set(CompTot,'String',sprintf('Total Compliance:%10.3f',-sum(CompE(:))/2))
set(CompWin,'String',sprintf('Window Compliance:%10.3f',-sum(CompE(use))/2))
%--------------------------------------------------------------------------

last={'Solid',num2str(iteration),0,MagSlide.Value};
lastW=W;    OldNodes=nodes;    p=0;
set(Xmax_Text,'String',W(2,1))
set(Xmin_Text,'String',W(1,1))
set(Ymax_Text,'String',W(2,2))
set(Ymin_Text,'String',W(1,2))
set(Zmax_Text,'String',W(2,3))
set(Zmin_Text,'String',W(1,3))

while(Done==0)
    pause(0.01)
    set(MagText,'String',sprintf('Mag:%3.1f',MagSlide.Value))
    if(~isequal(lastW,W))
        lastW=W;    p=1;
        set(Xmax_Text,'String',W(2,1))
        set(Xmin_Text,'String',W(1,1))
        set(Ymax_Text,'String',W(2,2))
        set(Ymin_Text,'String',W(1,2))
        set(Zmax_Text,'String',W(2,3))
        set(Zmin_Text,'String',W(1,3))
        %xlim(W(:,1));   ylim(W(:,2));   zlim(W(:,3));
        %Cent=nodes(elements(:,1),:)+Esize/2;
        use=find(Cent(:,1)>=W(1,1)&Cent(:,1)<=W(2,1)&...
            Cent(:,2)>=W(1,2)&Cent(:,2)<=W(2,2)&...
            Cent(:,3)>=W(1,3)&Cent(:,3)<=W(2,3));
        set(CompTot,'String',sprintf('Total Compliance:%10.3f',-sum(CompE(:))/2))
        set(CompWin,'String',sprintf('Window Compliance:%10.3f',-sum(CompE(use))/2))
        set(VolTot,'String',sprintf('Total Volume Fraction:%1.3f',prod(Esize)*sum(struc(map))/volTot))
        set(Volfrac,'String',sprintf('Window Volume Fraction:%1.3f',prod(Esize)*sum(struc(map(use)))/volTot))
        set(Intfrac,'String',sprintf('Window Interior Volume Fraction:%1.3f',prod(Esize)*sum(struc(map(setdiff(use,boundary))))/volTot))
    end
    if(~strcmp(itr.String,last{2}))
        p=1;
        iteration=max(1,min([str2num(itr.String),(numel(mat_files)-1)]));
        set(itr,'String',num2str(iteration));
        load([folder,'/Iteration',num2str(iteration)])
        OldNodes=nodes;
        Esize=max(nodes(elements(1,:),:))-min(nodes(elements(1,:),:));
        VonMises=zeros(size(elements,1),1);
        ind=find(struc(map));
        for(e=1:sum(struc(map)))
            stress=C*B*U(3*repelem(elements(ind(e),:),1,3)-repmat([2,1,0],1,8));
            VonMises(ind(e))=sqrt(sum((stress(1:3)-stress([2,3,1])).^2)+6*sum(stress(4:6).^2))/sqrt(2);
        end
        if(CalcComp==1)
            CompE=zeros(size(elements,1),1);
            dof=3*repelem(elements,1,3)-repmat([2,1,0],1,8);
            for(e=1:size(elements,1))
                CompE(e)=-max(struc(map(e)),0.0001)*U(dof(e,:))'*ke*U(dof(e,:));
            end
        end
        Cent=nodes(elements(:,1),:)+Esize/2;
        use=find(Cent(:,1)>=W(1,1)&Cent(:,1)<=W(2,1)&...
            Cent(:,2)>=W(1,2)&Cent(:,2)<=W(2,2)&...
            Cent(:,3)>=W(1,3)&Cent(:,3)<=W(2,3));
        set(CompTot,'String',sprintf('Total Compliance:%10.3f',-sum(CompE(:))/2))
        set(CompWin,'String',sprintf('Window Compliance:%10.3f',-sum(CompE(use))/2))
        set(VolTot,'String',sprintf('Total Volume Fraction:%1.3f',prod(Esize)*sum(struc(map))/volTot))
        set(Volfrac,'String',sprintf('Window Volume Fraction:%1.3f',prod(Esize)*sum(struc(map(use)))/volTot))
        set(Intfrac,'String',sprintf('Window Interior Volume Fraction:%1.3f',prod(Esize)*sum(struc(map(setdiff(use,boundary))))/volTot))
    end
    if(DeflecTB.Value~=last{3}||(MagSlide.Value~=last{4}&&DeflecTB.Value==1))
        p=1;
    end
    
    if(p==1||~strcmp(P.SelectedObject.String,last{1}))
        cla
        colorbar('off');    legend('off');
        nodes=OldNodes+DeflecTB.Value*MagSlide.Value*reshape(U,3,[])';
        if(strcmp(P.SelectedObject.String,'Solid'))
            plotstructure(elements(use,:),nodes,struc,map(use));
            lgd=legend('Solid');
            lgd.Position=[0.85,0.85,0.1,0.1];
        elseif(strcmp(P.SelectedObject.String,'Transparent'))
            plottrans(elements(use,:),nodes,struc,map(use),find(ismember(use,boundary)));
        elseif(strcmp(P.SelectedObject.String,'Stress'))
            plotstress(elements(use,:),nodes,struc,map(use),VonMises(use),Yield);
        else
            plotvoid(elements,nodes,struc);
            hold on;
            plotSTL('Rotated Irregular Pressure Vessel.STL');
            lgd=legend('Void');
            lgd.Position=[0.85,0.85,0.1,0.1];
        end
        last={P.SelectedObject.String,itr.String,DeflecTB.Value,MagSlide.Value};
        p=0;
    end
end
disp('Done')

%X Max Functions-----------------------------------------------------------
    function XmaxDPushed(scr,event)
        global Domain W Esize
        W(2,1)=min(Domain(2,1),max(W(1,1)+Esize(1),W(2,1)-Esize(1)));
    end

    function XmaxUPushed(scr,event)
        global Domain W Esize
        W(2,1)=min(Domain(2,1),max(W(1,1)+Esize(1),W(2,1)+Esize(1)));
    end
%--------------------------------------------------------------------------

%X Min Functions-----------------------------------------------------------
    function XminDPushed(scr,event)
        global Domain W Esize
        W(1,1)=max(Domain(1,1),min(W(2,1)-Esize(1),W(1,1)-Esize(1)));
    end

    function XminUPushed(scr,event)
        global Domain W Esize
        W(1,1)=max(Domain(1,1),min(W(2,1)-Esize(1),W(1,1)+Esize(1)));
    end
%--------------------------------------------------------------------------

%Y Max Functions-----------------------------------------------------------
    function YmaxDPushed(scr,event)
        global Domain W Esize
        W(2,2)=min(Domain(2,2),max(W(1,2)+Esize(2),W(2,2)-Esize(2)));
    end

    function YmaxUPushed(scr,event)
        global Domain W Esize
        W(2,2)=min(Domain(2,2),max(W(1,2)+Esize(2),W(2,2)+Esize(2)));
    end
%--------------------------------------------------------------------------

%Y Min Functions-----------------------------------------------------------
    function YminDPushed(scr,event)
        global Domain W Esize
        W(1,2)=max(Domain(1,2),min(W(2,2)-Esize(2),W(1,2)-Esize(2)));
    end

    function YminUPushed(scr,event)
        global Domain W Esize
        W(1,2)=max(Domain(1,2),min(W(2,2)-Esize(2),W(1,2)+Esize(2)));
    end
%--------------------------------------------------------------------------

%Z Max Functions-----------------------------------------------------------
    function ZmaxDPushed(scr,event)
        global Domain W Esize
        W(2,3)=min(Domain(2,3),max(W(1,3)+Esize(3),W(2,3)-Esize(3)));
    end

    function ZmaxUPushed(scr,event)
        global Domain W Esize
        W(2,3)=min(Domain(2,3),max(W(1,3)+Esize(3),W(2,3)+Esize(3)));
    end
%--------------------------------------------------------------------------

%Z Min Functions-----------------------------------------------------------
    function ZminDPushed(scr,event)
        global Domain W Esize
        W(1,3)=max(Domain(1,3),min(W(2,3)-Esize(3),W(1,3)-Esize(3)));
    end

    function ZminUPushed(scr,event)
        global Domain W Esize
        W(1,3)=max(Domain(1,3),min(W(2,3)-Esize(3),W(1,3)+Esize(3)));
    end
%--------------------------------------------------------------------------


%Plot type callbacks--------------------------------SolidPushed
    function SPushed(scr,event)
        %save figure
        num=numel(dir('*Cross*fig*'))+1;
        savefig(gcf,[pwd,'\','Cross-Section',num2str(num)]);
    end

    function DPushed(scr,event)
        global Done
        Done=1;
    end











