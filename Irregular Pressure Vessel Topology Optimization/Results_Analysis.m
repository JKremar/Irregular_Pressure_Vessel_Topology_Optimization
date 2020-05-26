clear all
close all
clc

folder='0.45finish23-May-2020 10,08,20';    %iIPV 30 v=0.45

C=reshape([0,1,0;1,2,1;0,1,0;1,2,1;2,3,2;1,2,1;0,1,0;1,2,1;0,1,0],3,3,3);

load([folder,'/Iteration0'])
mat_files=dir([folder,'/*.mat']);
max_itr=numel(mat_files)-1;
LSFsize=size(lsf);
fprintf('Analyzing %d iterations \n',max_itr);
Domain=[min(nodes);max(nodes)];
Esize=max(nodes(elements(1,:),:))-min(nodes(elements(1,:),:));
numelem=size(elements,1);
volTot=prod(Esize)*numelem;
InteriorLSF=setdiff(1:numel(lsf),bearing);



LSF_all=zeros([LSFsize,max_itr]);
Shape_all=zeros([LSFsize,max_itr]);
SensTot_all=zeros([LSFsize,max_itr]);
U_all=cell(max_itr,1);
ShapeError=zeros(4,max_itr);
SensError=zeros(4,max_itr);
LSF_data=zeros(3,max_itr);
Penalty_all=zeros(max_itr,1);
PenCalc=zeros(max_itr,1);
La_all=zeros(max_itr,1);
La2_all=zeros(max_itr,1);
Comp=zeros(max_itr,1);
Vol=zeros(max_itr,1);
Pen=zeros(max_itr,1);
dt=zeros(max_itr,1);
Meshsize=zeros(max_itr,1);
remesh=[];

for(i=1:max_itr)
    load([folder,'/Iteration',num2str(i)])
    Esize=max(nodes(elements(1,:),:))-min(nodes(elements(1,:),:));
    Meshsize(i)=Esize(1);
    if(Meshsize(max(1,i-1))~=Meshsize(i))
        remesh=[remesh,i-1];
    end
    
    %map for the iteration
    ind=round(nodes(elements(:,1),:)./Esize+1-Domain(1,:)./Esize);
    map=sub2ind(size(struc),ind(:,1),ind(:,2),ind(:,3));
    %----------------------
    %CompE for the iteration
    CompE=[];
    dof=3*repelem(elements,1,3)-repmat([2,1,0],1,8);
    for(e=1:size(elements,1))
        CompE(e)=-max(struc(map(e)),0.0001)*U(dof(e,:))'*ke*U(dof(e,:));
    end
    %-----------------------
    
    
    
    Shape_all(:,:,:,i)=shapeSens;
    SensTot_all(:,:,:,i)=SensTotal;
    Comp(i)=-sum(CompE(:));
    Vol(i)=prod(Esize)*sum(struc(map))/volTot;
    LSF_all(:,:,:,i)=lsf;
    U_all{i}=U;
    La_all(i)=La;
    La2_all(i)=La2;
    Penalty_all(i)=Penalty;
    PenCalc(i)=(La)*(Vol(i)-volReq);
    LSF_data(1,i)=mean(lsf(InteriorLSF));
    LSF_data(2,i)=mean(lsf(InteriorLSF))-min(lsf(InteriorLSF));
    LSF_data(3,i)=max(lsf(InteriorLSF))-mean(lsf(InteriorLSF));
    ShapeError(1,i)=mean(shapeSens(InteriorLSF));
    ShapeError(2,i)=std(shapeSens(InteriorLSF));
    ShapeError(3,i)=max(shapeSens(InteriorLSF));
    ShapeError(4,i)=min(shapeSens(InteriorLSF));
    P(i)=PID(1)*(Vol(i)-volReq);
    I(i)=PID(2)*((sum(Vol(max(1,i-4):i))/numel(max(1,i-4):i))-volReq);
    D(i)=PID(3)*(2*Vol(i)-Vol(max(1,i-1))-volReq);
    SensTotal=convn(padarray(SensTotal,[1,1,1],'replicate'),1/8*C,'valid');
    SensTotal(bearing)=0;
    dt(i)=Esize(1)*0.1/max(abs(SensTotal(:)));
    SensError(1,i)=mean(SensTotal(InteriorLSF));
    SensError(2,i)=std(SensTotal(InteriorLSF));
    SensError(3,i)=max(SensTotal(InteriorLSF));
    SensError(4,i)=min(SensTotal(InteriorLSF));
end
SidebySide=[Comp,Vol,Penalty_all,SensError(1,:)',dt,La_all];
x=1:i;
disp('Plotting')

mkdir(['Analysis',folder])
PIDplot=figure;
hold on
plot(x,P,x,I,x,D)
yyaxis right
plot(x,La_all,x,La2_all)
title('PID Terms')
legend({'Proportional Term','Integral Term','Derivative Term'});
xlabel('Iteration')
saveas(PIDplot,[pwd,'\Analysis',folder,'\','PID_Terms']);



ShapeE=figure;
hold on
errorbar([1:i],ShapeError(1,:),ShapeError(2,:))
plot(x,ShapeError(3,:))
xlabel('Iteration')
title('ShapeSens Error Plot')
legend({'Mean and STD','Maximum'})
saveas(ShapeE,[pwd,'\Analysis',folder,'\','ShapeError']);

SensE=figure;
hold on
errorbar([1:i],SensError(1,:),SensError(2,:))
plot(x,SensError(3,:),x,SensError(4,:))
xlabel('Iteration')
title('SensTotal Error Plot')
legend({'Mean and STD','Maximum'})
saveas(SensE,[pwd,'\Analysis',folder,'\','SensError']);


LSFdistribution=figure;
errorbar([1:i],LSF_data(1,:),LSF_data(2,:),LSF_data(3,:))
title('LSF Error Plot')
xlabel('Iteration')
legend({'Mean and Max/Min'})
saveas(LSFdistribution,[pwd,'\Analysis',folder,'\','LSFdistribution']);


Penalty=figure;
plot(Penalty_all)
hold on
plot(PenCalc)
title('Penalty Term')
xlabel('Iteration')
legend({'Used Penalty','Original Formulation'},'Location','southwest')
saveas(Penalty,[pwd,'\Analysis',folder,'\','Penalty']);
y1=get(gca,'ylim');



TimeStep=figure;
plot(dt)
title('dt')
xlabel('Iteration')
saveas(TimeStep,[pwd,'\Analysis',folder,'\','TimeStep']);


Compliance=figure;
plot(Comp)
set(gca,'YScale','log')
title('Compliance')
xlabel('Iteration')
ylim([min(mean(Comp)-std(Comp),Comp(end)),mean(Comp)+std(Comp)]);
saveas(Compliance,[pwd,'\Analysis',folder,'\','Compliance']);


Volume=figure;
plot(Vol)
hold on
xlabel('Iteration')
ylabel('Volume Fraction')
plot(volReq*ones(1,numel(Vol)))
title('Volume Fraction')
%xlim([45,60])
legend({'Iteration Volume','Goal Volume'},'Location','southwest')
saveas(Volume,[pwd,'\Analysis',folder,'\','Volume']);


ObjCon=figure;
subplot(2,1,1)
plot(Comp)
set(gca,'YScale','log')
title('Compliance')
ylim([min(mean(Comp)-std(Comp),Comp(end)),mean(Comp)+std(Comp)])
xlabel('Iteration')
subplot(2,1,2)
plot(Vol)
hold on
plot(volReq*ones(1,numel(Vol)))
title('Volume Fraction')
legend({'Iteration Volume','Goal Volume'},'Location','northeast')
xlabel('Iteration')
saveas(ObjCon,[pwd,'\Analysis',folder,'\','ComplianceAndVolume']);

summary=figure;
subplot(3,1,1)
hold on
plot(x,P,x,I,x,D)
yyaxis right
plot(x,La_all,x,La2_all)
title('PID Terms')
legend({'Proportional Term','Integral Term','Derivative Term'});
xlabel('Iteration')
subplot(3,1,2)
yyaxis left
plot(x,Vol,x,volReq*ones(1,numel(Vol)))
yyaxis right
plot(x,Comp)
set(gca,'YScale','log')
ylim([min(mean(Comp)-std(Comp),Comp(end)),mean(Comp)+std(Comp)])
title('Volume and Compliance')
legend({'Volume','Goal Volume','Compliance'});
xlabel('Iteration')
subplot(3,1,3)
plot(Penalty_all)
hold on
plot(PenCalc)
title('Penalty Term')
xlabel('Iteration')
legend({'Used Penalty','Original Formulation'},'Location','southwest')
saveas(summary,[pwd,'\Analysis',folder,'\','Summary Plot']);

disp('done')



