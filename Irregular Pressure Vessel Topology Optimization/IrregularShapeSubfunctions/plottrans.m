function [fig] = plottrans(elements,nodes,struc,map,boundary)


Eouter=elements(boundary,:);
Ecenter=elements(setdiff(find(struc(map)),boundary),:);
s(:,:,1) = Ecenter(:,[1,4,3,2]);
s(:,:,2) = Ecenter(:,[1,2,6,5]);
s(:,:,3) = Ecenter(:,[2,3,7,6]);
s(:,:,4) = Ecenter(:,[3,4,8,7]);
s(:,:,5) = Ecenter(:,[4,1,5,8]);
s(:,:,6) = Ecenter(:,[5,6,7,8]);
o(:,:,1) = Eouter(:,[1,4,3,2]);
o(:,:,2) = Eouter(:,[1,2,6,5]);
o(:,:,3) = Eouter(:,[2,3,7,6]);
o(:,:,4) = Eouter(:,[3,4,8,7]);
o(:,:,5) = Eouter(:,[4,1,5,8]);
o(:,:,6) = Eouter(:,[5,6,7,8]);


for(i=1:6)
    b=patch('Vertices',nodes,'Faces',o(:,:,i));
    set(b,'facecolor',[0.8485,0.49959,0.17446],'FaceAlpha',0.1,'edgecolor','none');
end
for(i=1:6)
    p=patch('Vertices',nodes,'Faces',s(:,:,i));
    set(p,'facecolor',[0.60551,0.38649,0.69569],'edgecolor','black');
end
lgd=legend([b,p],'Boundary','Solid')
lgd.Position=[0.85,0.85,0.1,0.1];


fig=gca;



end