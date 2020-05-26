function [fig] = plotstructure(elements,nodes,structure,map)

E=elements(find(structure(map)),:);

s(:,:,1) = E(:,[1,4,3,2]);
s(:,:,2) = E(:,[1,2,6,5]);
s(:,:,3) = E(:,[2,3,7,6]);
s(:,:,4) = E(:,[3,4,8,7]);
s(:,:,5) = E(:,[4,1,5,8]);
s(:,:,6) = E(:,[5,6,7,8]);

for(i=1:6)
    p=patch('Vertices',nodes,'Faces',s(:,:,i));
    set(p,'facecolor',[0.9290, 0.6940, 0.1250],'edgecolor','black','FaceLighting','gouraud','AmbientStrength',0.5);
end
camlight left; lighting phong;
fig=gca;

end
