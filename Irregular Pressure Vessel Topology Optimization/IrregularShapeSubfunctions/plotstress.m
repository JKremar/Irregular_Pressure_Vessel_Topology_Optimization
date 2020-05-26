function [fig] = plotstress(elements,nodes,structure,map,stress,yield)


E=elements(find(structure(map)),:);

s(:,:,1) = E(:,[1,4,3,2]);
s(:,:,2) = E(:,[1,2,6,5]);
s(:,:,3) = E(:,[2,3,7,6]);
s(:,:,4) = E(:,[3,4,8,7]);
s(:,:,5) = E(:,[4,1,5,8]);
s(:,:,6) = E(:,[5,6,7,8]);

colormap(jet)
caxis([0,yield])
C=min(yield,stress(find(structure(map))));%(find(structure(map)))/yield);
fprintf('Max VonMises Stress: %10.2f Average Stress: %10.2f Yield: %10.2f\n',max(stress),mean(stress),yield);

for(i=1:6)
    p=patch('Vertices',nodes,'Faces',s(:,:,i),'FaceVertexCData',C,'FaceColor','flat');
    set(p,'FaceLighting','gouraud','AmbientStrength',0.5);
end
camlight left; lighting phong;
cb=colorbar('Position',[0.95,0.25,0.025,0.65],'AxisLocation','in');
cb.Ticks=linspace(0,yield,6);
cb.TickLabels=strsplit([num2str(linspace(0,yield,6)),'\newline{Yield}']);

fig=gca;

end