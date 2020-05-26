function [fig] = plotSTL(file)

addpath([pwd,'\MakeMeshSubfunctions'])

[faces,vertices] = readSTL(file,'inches');

vertices=reshape(permute(vertices,[2,1,3]),3,[])';
vertices=vertices-min(vertices);

faces=reshape(1:size(vertices,1)/3,3,[])';
p=patch('Vertices',vertices,'Faces',faces);
set(p,'facecolor',[0.60551,0.38649,0.69569],'FaceAlpha',0.1,'edgecolor','black');


fig=gca;



end
