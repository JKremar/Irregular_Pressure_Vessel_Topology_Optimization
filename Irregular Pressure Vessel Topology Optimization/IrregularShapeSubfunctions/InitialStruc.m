function [struc,Esize,map,noF,exterior] =InitialStruc(elements,nodes,boundary,init)
%Evaluated initial values for optimization
%

Domain=[min(nodes);max(nodes)];
Esize=max(nodes(elements(1,:),:))-min(nodes(elements(1,:),:));
StrucSize=round((Domain(2,:)-Domain(1,:))./Esize);
struc=ones(StrucSize);
ind=round(nodes(elements(:,1),:)./Esize+1-Domain(1,:)./Esize);
map=sub2ind(StrucSize,ind(:,1),ind(:,2),ind(:,3));  %map is a list for each element, which struc index is used
void=zeros(init(1,1)*init(3,1)+init(2,1)*(init(3,1)-1),...
    init(1,2)*init(3,2)+init(2,2)*(init(3,2)-1),...
    init(1,3)*init(3,3)+init(2,3)*(init(3,3)-1));
void([0:init(3,1)-1]'*(init(1,1)+init(2,1))+[1:init(1,1)],...
    [0:init(3,2)-1]'*(init(1,2)+init(2,2))+[1:init(1,2)],...
    [0:init(3,3)-1]'*(init(1,3)+init(2,3))+[1:init(1,3)])=1;
vs=size(void);
bounds=round(mean(nodes)./Esize-vs/2);
void=void(max(1,2-bounds(1)):min(vs(1),StrucSize(1)-bounds(1)-1),...
    max(1,2-bounds(2)):min(vs(2),StrucSize(2)-bounds(2)-1),...
    max(1,2-bounds(3)):min(vs(3),StrucSize(3)-bounds(3)-1));
bounds=[max(2,bounds);max(2,bounds)+size(void)-1];
struc(bounds(1,1):bounds(2,1),bounds(1,2):bounds(2,2),bounds(1,3):bounds(2,3))=...
    max(0,struc(bounds(1,1):bounds(2,1),bounds(1,2):bounds(2,2),bounds(1,3):bounds(2,3))-void);
struc(map(boundary))=1;
struc(setdiff(1:prod(StrucSize),map))=1;


bn=elements(boundary,:);
bn=unique(bn(:));
noFnodes=bn(find(sum(bn'==elements(:))<8));
noF=reshape(3*noFnodes'-[2;1;0],1,[]);
exterior=(setdiff(1:numel(struc),map))';

end

