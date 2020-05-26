function [fig] = plotvoid(elements,nodes,structure)

Esize=max(nodes(elements(1,:),:))-min(nodes(elements(1,:),:));
void=find(~structure);
[r,c,p]=ind2sub(size(structure),void);
cent=[r,c,p].*Esize-0.5*Esize;
Vnodes=permute(cent,[3,2,1])+0.5.*Esize.*[-1,-1,-1;1,-1,-1;1,1,-1;-1,1,-1;-1,-1,1;1,-1,1;1,1,1;-1,1,1];
Vnodes=reshape(permute(Vnodes,[2,1,3]),3,[])';
E=reshape(1:size(Vnodes,1),8,[])';

s(:,:,1) = E(:,[1,4,3,2]);
s(:,:,2) = E(:,[1,2,6,5]);
s(:,:,3) = E(:,[2,3,7,6]);
s(:,:,4) = E(:,[3,4,8,7]);
s(:,:,5) = E(:,[4,1,5,8]);
s(:,:,6) = E(:,[5,6,7,8]);


for(i=1:6)
    p=patch('Vertices',Vnodes,'Faces',s(:,:,i));
    set(p,'facecolor','yellow','edgecolor','black');
end

fig=gca;



end
