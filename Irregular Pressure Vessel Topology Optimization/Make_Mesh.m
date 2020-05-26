close all
clear all
clc
addpath([pwd,'\MakeMeshSubfunctions'])

voxelsize=[0.0625,0.0625,0.0625];       %element pixelsize in the x, y, and z direction 

if(numel(voxelsize)==1)
    voxelsize(1:3)=voxelsize;
end


%[faces,vertices] = readSTL('Rotated Irregular Pressure Vessel.STL','inches');
[faces,vertices] = readSTL('Oxygen Tank Hallow.STL','inches');

ranges=[min(min(vertices),[],3);max(max(vertices),[],3)];
x_centroids=ranges(1,1)+0.5*voxelsize(1):voxelsize(1):ranges(2,1);
y_centroids=ranges(1,2)+0.5*voxelsize(2):voxelsize(2):ranges(2,2);
z_centroids=ranges(1,3)+0.5*voxelsize(3):voxelsize(3):ranges(2,3);
xyfaces=find(faces(:,3)~=0);
disp('STL file read')

[X,Y]=meshgrid(x_centroids,y_centroids);
numxyzf=[fliplr(size(X)),numel(z_centroids),numel(xyfaces)];
cells=-1*ones(numxyzf([2,1,3]));
for(f=1:numxyzf(4))
    p=vertices(:,1:2,xyfaces(f));
    intersects=inpolygon(X,Y,p(:,1),p(:,2));
    A=0.5*det([[1;1;1],p]);
    Plane=zeros(2);
    for(i=1:3)
        up=rem(i,3)+1;
        down=3-rem(4-i,3);
        Plane=Plane+0.5*vertices(i,3,xyfaces(f))*[0,p(down,1)-p(up,1);p(up,2)-p(down,2),p(up,1)*p(down,2)-p(down,1)*p(up,2)]/A;
    end
    I=find(intersects);
    if(~isempty(I))
        z_int=poly2Deval(Plane,[X(I),Y(I)]);
        mult=-1*(z_centroids>=z_int)+(z_centroids<z_int);
        %[z_int,mult]
        t=0;
        for(i=1:numel(I))
            %mult=-1*(z_centroids>=z_int(i))+(z_centroids<z_int(i))
            [r,c]=ind2sub([numxyzf(2),numxyzf(3)],I(i));
            cells(r,c,:)=cells(r,c,:).*permute(mult(i,:),[1,3,2]);
            %cells(r,c,:)=cells(r,c,:).*permute(mult,[1,3,2]);
        end
    end
    if(~mod(f,250))
        fprintf('evaluated %d of %d faces \n',f,numxyzf(4));
    end
end



disp('generating mesh')
cells=permute(cells,[2,1,3]);
%outer=(padarray(cells,[1,1,1],-1)==1).*(bwdist(~(padarray(cells,[1,1,1],-1)/2+0.5))<2);
%boundary=

outer2=(cells==1).*(convn(cells,ones(3,3,3),'same')<27);
outer2(cells(:)==-1)=[];
boundary=nonzeros(outer2(:)'.*(1:nnz(cells==1)));


nelx=numxyzf(1);    nely=numxyzf(2);    nelz=numxyzf(3);
Elements=1:nelx*nely*nelz;
n1z=floor((Elements-1)/(nelx*nely));
n1x=rem((Elements-(nelx*nely).*floor((Elements-1)/(nelx*nely))-1),nelx);
n1y=floor((Elements-(nelx*nely).*floor((Elements-1)/(nelx*nely))-1)/nelx);
Relative=[0;1;nelx+2;nelx+1;...
    (nelx+1)*(nely+1);(nelx+1)*(nely+1)+1;(nelx+1)*(nely+1)+nelx+2;(nelx+1)*(nely+1)+nelx+1];
Elements=(1+n1x+n1y*(nelx+1)+n1z*(nelx+1)*(nely+1))'+Relative';
[Nodes(:,1),Nodes(:,2),Nodes(:,3)]=ind2sub([nelx+1,nely+1,nelz+1],1:(nelx+1)*(nely+1)*(nelz+1));
Nodes=voxelsize.*(Nodes-[1,1,1]);


Elements=(cells(:)==1).*Elements;   %find on elements
Elements((Elements(:,1)==0),:)=[];  %remove off elements

elements=zeros(size(Elements));
nodes=zeros(size(Nodes));
N=1;
while(sum(elements(:)==0)>0)
    [c,r]=find(elements'==0,1);
    ind=find(Elements==Elements(r,c));
    elements(ind)=N;
    nodes(N,:)=Nodes(Elements(r,c),:);
    N=N+1;
    if(~mod(N,1000))
        fprintf('meshing node %d of %d \n',sum(elements(:)~=0),numel(Elements));
    end
end
nodes(N:end,:)=[];

if(1==2)
    disp('Finding Constraints')
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
else
    fix=[];
end

disp('plotting')
meshplot(elements, nodes, boundary);    axis equal
view(44,-18)
xlabel('x axis')
ylabel('y axis')
zlabel('z axis')
%meshplotlayer(elements, nodes, Boundary);
disp('done')

%to safe information----------------------------------------
%save('filename','variables',...)
save('SphericalTank0625','elements','nodes','boundary','fix')
%-----------------------------------------------------------

