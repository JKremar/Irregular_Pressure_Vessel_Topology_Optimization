function [struc,elements,nodes,map,boundary,noF,sX,sY,sZ,exterior] = remesh(lsf,LSFcoord,Esize,file)
%Remeshes material domain excluding void regions

[faces,vertices]=readSTL(file,'inches');

ranges=[min(min(vertices),[],3);max(max(vertices),[],3)];
x_centroids=ranges(1,1)+0.5*Esize(1):Esize(1):ranges(2,1);
y_centroids=ranges(1,2)+0.5*Esize(2):Esize(2):ranges(2,2);
z_centroids=ranges(1,3)+0.5*Esize(3):Esize(3):ranges(2,3);
xyfaces=find(faces(:,3)~=0);

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
        for(i=1:numel(I))
            [r,c]=ind2sub([numxyzf(2),numxyzf(3)],I(i));
            cells(r,c,:)=cells(r,c,:).*permute(mult(i,:),[1,3,2]);
        end
    end
end
cells=permute(cells,[2,1,3]);
exterior=find(cells==-1);

Domain=ranges-ranges(1,:);
[sX,sY,sZ]=meshgrid(Esize(1)/2:Esize(1):Domain(2,1),...
    Esize(2)/2:Esize(2):Domain(2,2),Esize(3)/2:Esize(3):Domain(2,3));
sX=permute(sX,[2,1,3]); sY=permute(sY,[2,1,3]); sZ=permute(sZ,[2,1,3]);
struc=griddata(LSFcoord(:,1),LSFcoord(:,2),LSFcoord(:,3),lsf,sX,sY,sZ)<=0;


Fullcells=padarray(cells,[1,1,1],-1);
Nf_1458=(Fullcells-circshift(Fullcells,[1,0,0]))==2;
Nf_2367=(Fullcells-circshift(Fullcells,[-1,0,0]))==2;
Nf_1256=(Fullcells-circshift(Fullcells,[0,1,0]))==2;
Nf_3478=(Fullcells-circshift(Fullcells,[0,-1,0]))==2;
Nf_1234=(Fullcells-circshift(Fullcells,[0,0,1]))==2;
Nf_5678=(Fullcells-circshift(Fullcells,[0,0,-1]))==2;


outer=(cells==1).*(convn(cells,ones(3,3,3),'same')<27);
struc(find(outer))=1;


nelx=numxyzf(1);    nely=numxyzf(2);    nelz=numxyzf(3);
Elements=1:nelx*nely*nelz;
n1z=floor((Elements-1)/(nelx*nely));
n1x=rem((Elements-(nelx*nely).*floor((Elements-1)/(nelx*nely))-1),nelx);
n1y=floor((Elements-(nelx*nely).*floor((Elements-1)/(nelx*nely))-1)/nelx);
Relative=[0;1;nelx+2;nelx+1;...
    (nelx+1)*(nely+1);(nelx+1)*(nely+1)+1;(nelx+1)*(nely+1)+nelx+2;(nelx+1)*(nely+1)+nelx+1];
Elements=(1+n1x+n1y*(nelx+1)+n1z*(nelx+1)*(nely+1))'+Relative';
[Nodes(:,1),Nodes(:,2),Nodes(:,3)]=ind2sub([nelx+1,nely+1,nelz+1],1:(nelx+1)*(nely+1)*(nelz+1));
Nodes=Esize.*(Nodes-[1,1,1]);


Elements=(cells(:)==1 & struc(:)==1).*Elements;   %find on elements
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
    if(~mod(N,5000))
        fprintf('meshing node %d of %d \n',sum(elements(:)~=0),numel(Elements));
    end
end
nodes(N:end,:)=[];

ind=round(nodes(elements(:,1),:)./Esize+1-Domain(1,:)./Esize);
map=sub2ind(numxyzf([1,2,3]),ind(:,1),ind(:,2),ind(:,3));  %map is a list for each element, which struc index is used
mapFull=sub2ind(numxyzf([1,2,3])+2,ind(:,1)+1,ind(:,2)+1,ind(:,3)+1);

boundary=find(outer(map));
noFnodes=[elements(find(Nf_1458(mapFull)),[1,4,5,8]);...
    elements(find(Nf_2367(mapFull)),[2,3,6,7]);...
    elements(find(Nf_1256(mapFull)),[1,2,5,6]);...
    elements(find(Nf_3478(mapFull)),[3,4,7,8]);...
    elements(find(Nf_1234(mapFull)),[1,2,3,4]);...
    elements(find(Nf_5678(mapFull)),[5,6,7,8])];
noF=reshape(3*unique(noFnodes(:)')-[2;1;0],1,[]);


end

