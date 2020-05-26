close all
clear all
clc
addpath([pwd,'\IrregularShapeSubfunctions'])
addpath([pwd,'\MakeMeshSubfunctions'])
%Attenpt to implement level-set topology optimization on a 3D structure with a
%pressure load being applied from a void in the center
%Calculates forces as outward normal for all void elements
disp('running...')
Title='RemeshIPVVol0.45';
%Material Parameters and Working Pressure------------
%Inconel718
E=29.5*10^6;    %psi
nu=0.29;
Yield=150*10^3; %psi
Pressure=3000;  %PSI
%----------------------------------------------------

%Geometry and Loading---------------
load('RotatedIPVmesh25.mat')   %imports 'elements' 'nodes' and 'boundary' from saved mesh file
load('RemeshStart25.mat')%load values of oldStruct,OldK,OldF to compare to
%-------------------------------

% Establish Level-Set parameters-----------------------
volReq=0.25; 
stepLength=2;
numReinit=3;
topWeight=0;
max_itr=200;
LSFspacing=0.375;
init=[3,3,3;2,2,2;10,10,10];    %edge length of initial void; gap between; repeated     
maxNodes=75000;
La=1/2;    La2=1;  alpha=1/0.95;
PID=[0.5,0.2,1];    relax=0;
%-------------------------------------------------------

%Initialization--------------------------------------------------------
i=1;    flag=0; mesh=1; band=0.15;  s=1;
%Initialize Struc----------------------------
Domain=[min(nodes);max(nodes)];
[struc,Esize,map,noF,exterior]=InitialStruc(elements,nodes,boundary,init);   %map is a list for each element, which struc index is used
[sX,sY,sZ]=meshgrid(Esize(1)/2:Esize(1):Domain(2,1),...
    Esize(2)/2:Esize(2):Domain(2,2),Esize(3)/2:Esize(3):Domain(2,3));
sX=permute(sX,[2,1,3]); sY=permute(sY,[2,1,3]); sZ=permute(sZ,[2,1,3]);
numelem=size(elements,1);   numnodes=size(nodes,1);
CompE=zeros(numelem,1);
volTot=prod(Esize)*numelem;
%Initialize LSF----------------------------
LSFsize=ceil(Domain(2,:)/LSFspacing)+1;
cent=mean(Domain);
lsfX=LSFspacing*(LSFsize(1)-1)*linspace(-0.5,0.5,LSFsize(1))+cent(1);
lsfY=LSFspacing*(LSFsize(2)-1)*linspace(-0.5,0.5,LSFsize(2))+cent(2);
lsfZ=LSFspacing*(LSFsize(3)-1)*linspace(-0.5,0.5,LSFsize(3))+cent(3);
[lsfX,lsfY,lsfZ]=meshgrid(lsfX,lsfY,lsfZ);
lsfX=permute(lsfX,[2,1,3]); lsfY=permute(lsfY,[2,1,3]); lsfZ=permute(lsfZ,[2,1,3]);
sdf=(~struc).*bwdist(struc)-struc.*bwdist(struc-1); %reinitialize LSF
lsf=griddata(sX,sY,sZ,double(sdf),lsfX,lsfY,lsfZ);
Nanind=find(isnan(lsf));
LSF2EleDist=(nodes(elements(:,1),1)'+Esize(1)/2'-lsfX(:)).^2+...
            (nodes(elements(:,1),2)'+Esize(2)/2'-lsfY(:)).^2+...
            (nodes(elements(:,1),3)'+Esize(3)/2'-lsfZ(:)).^2;
%LSF2StrucDist=(sX(:)'-lsfX(:)).^2+(sY(:)'-lsfY(:)).^2+(sZ(:)'-lsfZ(:)).^2;
[d,id]=min(LSF2EleDist(Nanind,:),[],2);
lsf(Nanind)=sdf(map(id))-d./Esize(1);
struc=griddata(lsfX,lsfY,lsfZ,lsf,sX,sY,sZ)<=0;
%Filter and Update Prep---------------------------------------
R=1.25*LSFspacing;
Hij=max(R-LSF2EleDist,0);

inside=setdiff(1:numelem,boundary);
bearing=find(sum(lsfX(:)'>=nodes(elements(inside,1),1) &...
    lsfY(:)'>=nodes(elements(inside,1),2) &...
    lsfZ(:)'>=nodes(elements(inside,1),3) &...
    lsfX(:)'<=nodes(elements(inside,7),1) &...
    lsfY(:)'<=nodes(elements(inside,7),2) &...
    lsfZ(:)'<=nodes(elements(inside,7),3))==0);
%-------------------------------------------------------------------------------------

%Loading and Boundary Conditions-------------------------------------------
Po=circshift(Esize,1).*circshift(Esize,-1)*Pressure/4;
dof=3*repelem(elements,1,3)-repmat([2,1,0],1,8);
[ke,B,C]=stiff3D(E,nu,Esize);
if(~isequal(struc,oldstruc))
    oldstruc=[];  oldK=[];  oldF=[];
end
if(~exist('fix'))
    Nborder=nodes(elements(boundary,:),:);
    D_best=0;
    for(b=1:size(Nborder,1))    %finds border points that would be best for coordinate oriented B.C.s
        N1=Nborder(b,:);
        Nx=Nborder(Nborder(:,2)==N1(2) & Nborder(:,3)==N1(3),:);
        Ny=Nborder(Nborder(:,1)==N1(1) & Nborder(:,3)==N1(3),:);
        Nz=Nborder(Nborder(:,1)==N1(1) & Nborder(:,2)==N1(2),:);
        [dx,ix]=max(abs(sum(N1-Nx,2)));
        [dy,iy]=max(abs(sum(N1-Ny,2)));
        [dz,iz]=max(abs(sum(N1-Nz,2)));
        if(dx^2+dy^2+dz^2>D_best)
            D_best=dx^2+dy^2+dz^2;
            fix=[N1;Nx(ix,:);Ny(iy,:);Nz(iz,:)];
        end
    end
end
[tf,ind]=ismember(fix,nodes,'rows');
if(sum(tf)~=4)
    disp('constraint error')
end
fixeddofs=nonzeros(reshape((3*ind-[2,1,0]).*[1,1,1;~eye(3)],[],1));
%--------------------------------------------------------------------------

%Save Initial-------------------------------------------------------------- 
folder=strcat(Title,strrep(datestr(datetime),':',','));
mkdir(folder);
save([pwd,'\',folder,'\','Iteration0'],'lsf','struc','La','alpha',...
    'init','volReq','ke','bearing','elements','nodes','boundary','map',...
    'max_itr','numReinit','Po','stepLength','PID','volTot','lsfX','lsfY','lsfZ')
clear LSF2EleDist;
disp(['Starting ',Title])
%--------------------------------------------------------------------------

while(flag==0)
    [U,oldK,oldF]=FEA_3DP6(struc,elements,map,ke,Po,noF,fixeddofs,oldstruc,oldK,oldF);
    %evaluate sensitivities of each element--------------------------------
    for(e=1:numelem)
        CompE(e)=-max(struc(map(e)),0.0001)*U(dof(e,:))'*ke*U(dof(e,:));
    end
    
    %Post Processing and Plotting------------------------------------------
    obj(i)=-sum(CompE(:));  
    vol(i)=prod(Esize)*sum(struc(map))/volTot;
    disp(['It.:' num2str(i) ' Compl.:' sprintf('%10.4f',obj(i)/2) ' Vol.: ' sprintf('%6.3f',vol(i))...
        '  La:' sprintf('%10.3f',La) '  LaPID:' sprintf('%10.5f',La2)])
    
    %check for convergence-------------------------------------------------
    if(i>5)
        if((abs(vol(i)-volReq)<0.005) && all(abs(obj(end)-obj(end-5:end-1))<0.03*abs(obj(end))))
            flag=1;
        end
        if(i>=max_itr)
            flag=2;
        end
    end
    
    %Update Procedure------------------------------------------------------     
    if(relax==0 && abs(vol(i)-volReq)<=0.035)                          %(max(abs(vol(i-4:i)-volReq))<0.05 && relax==0)
        relax=1;    %Stop relaxed penalty if within volume band (0.15)
        s=0.3;
        shapeSens=reshape((Hij*CompE)./max(sum(Hij,2),0.0001),LSFsize);
        SensTotal=(shapeSens/max(abs(shapeSens(:))));
        Con=-0.75:0.005:0.75;
        V=zeros(numel(Con),1);
        for(c=1:numel(Con))
            [newlsf]=updatestep3(lsf,SensTotal+Con(c),stepLength,bearing,Esize(1));
            newstruc=(griddata(lsfX,lsfY,lsfZ,newlsf,sX,sY,sZ)<=0);
            newstruc(map(boundary))=1; newstruc(exterior)=1;
            add=setdiff(find((newstruc-struc)==1),[map;exterior]);
            newmap=[map;add];
            V(c)=prod(Esize)*(sum(newstruc(newmap)))/volTot;
        end
        Control(i-1)=Con(find(V>=vol(i),1,'last'));
    end
    if(relax==0)    %Execute relaxed penalty
        La=alpha*La;
        Penalty=La*(vol(i)-volReq);
%         Control=[];
%         Control(i)=Penalty;
    else
        if(max(vol(max(1,i-5):i))-min(vol(max(1,i-5):i))<0.005 && i>5)
            La2=(alpha^2)*La2;  %Update Lagrange multiplier on PID if volume hasn't changed
        end
        Control(i)=La2*PID*[(vol(i)-volReq);...
            ((sum(vol(max(1,i-4):i))/numel(max(1,i-4):i))-volReq);...
            (2*vol(i)-vol(max(1,i-1))-volReq)];
        Penalty=sum(Control);
    end
    shapeSens=reshape((Hij*CompE)./max(sum(Hij,2),0.0001),LSFsize);
    SensTotal=(shapeSens/max(abs(shapeSens(:))))+Penalty;
    
    %Save values every iteration-------------------------------------------
    save([pwd,'\',folder,'\','Iteration',num2str(i)],'lsf','struc','U',...
        'La','La2','shapeSens','SensTotal','Penalty','oldstruc',...
        'oldK','oldF','nodes','elements','map','CompE','exterior','boundary')
    %----------------------------------------------------------------------
    oldlsf=lsf;
    [lsf]=updatestep3(lsf,SensTotal,stepLength,bearing,s*Esize(1));
    oldstruc=struc;
    struc=(griddata(lsfX,lsfY,lsfZ,lsf,sX,sY,sZ)<=0);
    struc(map(boundary))=1; struc(exterior)=1;
    %----------------------------------------------------------------------
    
    add=setdiff(find((struc-oldstruc)==1),[map;exterior]);
    newmap=[map;add];
    if((prod(Esize)*(sum(struc(newmap)))/volTot)<(volReq-0.04))
        disp('Stepped Back')
        lsf=oldlsf;
        struc=(griddata(lsfX,lsfY,lsfZ,lsf,sX,sY,sZ)<=0);
        struc(map(boundary))=1; struc(exterior)=1;
    end
            
    %Prep next iteration---------------------------------------------------
    if(mesh>=5)
        if(mesh>=8 && max(abs(vol(i-4:i)-volReq))<band && Esize(1)>(3/32))
            mesh=0;   disp('option1');
            band=0.8*band;
            Esize=max(0.75*Esize,1/8);
        elseif(sum(struc(map)==0)>0.5*numelem || numnodes>150000)
            mesh=0;   disp('option2');
            band=0.15;
            Esize=repelem((prod(Esize)*(sum(struc(map)+...
                numel(setdiff(find((struc-oldstruc)==1),[map;exterior]))))/(1.2*numelem))^(1/3),3);
        end
    end
    
    
    if(mesh==0) %Remesh
        fprintf('remeshing with element size: %2.5f\n',Esize(1));
        [struc,elements,nodes,map,boundary,noF,sX,sY,sZ,exterior]=remesh(lsf,...
            [lsfX(:),lsfY(:),lsfZ(:)],Esize,'Rotated Irregular Pressure Vessel.STL');
        numelem=size(elements,1);   numnodes=size(nodes,1);
        fprintf('meshing complete with %d elements and %d nodes\n',numelem,numnodes);
        while(numnodes>160000)   %retry if too many nodes
            Esize=Esize*1.05;
            fprintf('remeshing with element size: %2.5f\n',Esize(1));
            [struc,elements,nodes,map,boundary,noF,sX,sY,sZ,exterior]=...
                remesh(lsf,[lsfX(:),lsfY(:),lsfZ(:)],Esize,'Rotated Irregular Pressure Vessel.STL');
            numelem=size(elements,1);   numnodes=size(nodes,1);
            fprintf('meshing complete with %d elements and %d nodes\n',numelem,numnodes);
        end
        LSF2EleDist=(nodes(elements(:,1),1)'+Esize(1)/2'-lsfX(:)).^2+...
            (nodes(elements(:,1),2)'+Esize(2)/2'-lsfY(:)).^2+...
            (nodes(elements(:,1),3)'+Esize(3)/2'-lsfZ(:)).^2;
        [d,id]=min(LSF2EleDist(Nanind,:),[],2);
        Hij=max(R-LSF2EleDist,0);
        Po=circshift(Esize,1).*circshift(Esize,-1)*Pressure/4;
        dof=3*repelem(elements,1,3)-repmat([2,1,0],1,8);
        [ke,B,C]=stiff3D(E,nu,Esize);
        oldstruc=[];  oldK=[];  oldF=[];
        [~,ind]=min((fix(:,1)-nodes(:,1)').^2+(fix(:,2)-nodes(:,2)').^2+(fix(:,3)-nodes(:,3)').^2,[],2);
        fixeddofs=nonzeros(reshape((3*ind-[2,1,0]).*[1,1,1;~eye(3)],[],1)); 
        clear LSF2EleDist;
        mesh=1;
    elseif((prod(Esize)*sum(struc(map))/volTot)>0.98) %If entire domain becomes solid revert back to initial configuration
        disp('Domain solid reverting to original discritization')
        La=(alpha)^5*La;    %Take a large step in La
        load('RotatedIPVmesh25.mat')
        load('RemeshStart25.mat')
        [struc,Esize,map,noF]=InitialStruc(elements,nodes,boundary,init);
        sdf=(~struc).*bwdist(struc)-struc.*bwdist(struc-1); %reinitialize LSF
        lsf=griddata(sX,sY,sZ,double(sdf),lsfX,lsfY,lsfZ);
        Nanind=find(isnan(lsf));
        LSF2EleDist=(nodes(elements(:,1),1)'+Esize(1)/2'-lsfX(:)).^2+...
            (nodes(elements(:,1),2)'+Esize(2)/2'-lsfY(:)).^2+...
            (nodes(elements(:,1),3)'+Esize(3)/2'-lsfZ(:)).^2;
        [d,id]=min(LSF2EleDist(Nanind,:),[],2);
        lsf(Nanind)=sdf(map(id))-d./Esize(1);
        struc=griddata(lsfX,lsfY,lsfZ,lsf,sX,sY,sZ)<=0;
        Hij=max(R-LSF2EleDist,0);
        Po=circshift(Esize,1).*circshift(Esize,-1)*Pressure/4;
        dof=3*repelem(elements,1,3)-repmat([2,1,0],1,8);
        [ke,B,C]=stiff3D(E,nu,Esize);
        if(~isequal(struc,oldstruc))
            oldstruc=[];  oldK=[];  oldF=[];
        end
        clear LSF2EleDist;
        mesh=1;
    else
        if(~mod(i,numReinit))   %reinitialize LSF
            sdf=(~struc).*bwdist(struc)-struc.*bwdist(struc-1); %reinitialize LSF
            lsf=griddata(sX,sY,sZ,double(sdf),lsfX,lsfY,lsfZ);
            lsf(Nanind)=sdf(map(id))-d./Esize(1);
            struc=griddata(lsfX,lsfY,lsfZ,lsf,sX,sY,sZ)<=0;
            struc(map(boundary))=1; struc(exterior)=1;
            clear sdf
        end
        mesh=mesh+1;    
        %Add elements if needed
        add=setdiff(find((struc-oldstruc)==1),[map;exterior]);
        if(~isempty(add))
            oldnumN=numnodes;
            mult=[-1,-1,-1;1,-1,-1;1,1,-1;-1,1,-1;-1,-1,1;1,-1,1;1,1,1;-1,1,1];
            Nfull=permute([sX(add),sY(add),sZ(add)],[3,2,1])+(Esize/2).*mult;
            Nall=reshape(permute(Nfull,[1,3,2]),[],3);
            [~,Nnum]=ismembertol(Nall,nodes,0.01*Esize(1),'ByRows',true);
            nodes=[nodes;Nall(Nnum==0,:)];
            Nnum(Nnum==0)=(numnodes+1):(numnodes+sum(Nnum==0));
            elements=[elements;reshape(Nnum,8,[])'];
            numnodes=size(nodes,1);
            numelem=size(elements,1);
            map=[map;add];
            dof=[dof;3*repelem(elements(end-numel(add)+1:end,:),1,3)-repmat([2,1,0],1,8)];
            L2ED=(nodes(elements(end-numel(add)+1:end,1),1)'+Esize(1)/2'-lsfX(:)).^2+...
                (nodes(elements(end-numel(add)+1:end,1),2)'+Esize(2)/2'-lsfX(:)).^2+...
                (nodes(elements(end-numel(add)+1:end,1),3)'+Esize(3)/2'-lsfX(:)).^2;
            Hij=[Hij,max(R-L2ED,0)];
            K=oldK; F=oldF;
            oldF=[oldF;zeros(3*(numnodes-oldnumN),1)];
            oldK=sparse(3*numnodes,3*numnodes);
            oldK(1:3*oldnumN,1:3*oldnumN)=K;
            for(a=1:numel(add))
                oldK(dof(end+1-a,:),dof(end+1-a,:))=oldK(dof(end+1-a,:),dof(end+1-a,:))+0.0001*ke;
            end
            fprintf('Added %d elements, new node total:%d\n',numel(add),numnodes)
            clear K a N Nnum add L2ED
        end
    end
    CompE=zeros(numelem,1);
    i=i+1;
    %----------------------------------------------------------------------
end
%End of Optimization

%Show Final Values
disp('done')


