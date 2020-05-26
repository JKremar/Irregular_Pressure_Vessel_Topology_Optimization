function [faces,vertices] = readSTL(filename,units)
%Reads an STL file found under the presribed filename and then filters and
%outputs the face normals and vertices of each face


fid=fopen(filename,'r');

Title=fgetl(fid);
f=0;
while(feof(fid)==0)
    line=fgetl(fid);
    if(contains(line,'facet normal'))
        f=f+1;
        v=0;
        faces(f,:)=str2num(line(17:end));
    elseif(contains(line,'vertex'))
        v=v+1;
        vertices(v,:,f)=str2num(line(17:end));
    end 
end
fclose(fid);


if(units=='inches')
    vertices=vertices/25.4;
end


end

