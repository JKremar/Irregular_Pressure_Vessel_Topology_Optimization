function [] = meshplot(elements, nodes, outer)
%Plots a 3D mesh


if(nargin<3)
    s(:,:,1) = elements(:,[1,4,3,2]);
    s(:,:,2) = elements(:,[1,2,6,5]);
    s(:,:,3) = elements(:,[2,3,7,6]);
    s(:,:,4) = elements(:,[3,4,8,7]);
    s(:,:,5) = elements(:,[4,1,5,8]);
    s(:,:,6) = elements(:,[5,6,7,8]);


    for i = 1:6
        p=patch('Vertices',nodes,'Faces',s(:,:,i));
        set(p,'facecolor','yellow','edgecolor','black');
    end
    
else
    boundary=elements(outer,:);
    elements=elements(setdiff(1:size(elements,1),outer),:);
    %nonzeros(outer(:).*elements);
    %elements=nonzeros(elements-(outer(:).*elements));
    s(:,:,1) = elements(:,[1,4,3,2]);
    s(:,:,2) = elements(:,[1,2,6,5]);
    s(:,:,3) = elements(:,[2,3,7,6]);
    s(:,:,4) = elements(:,[3,4,8,7]);
    s(:,:,5) = elements(:,[4,1,5,8]);
    s(:,:,6) = elements(:,[5,6,7,8]);
    o(:,:,1) = boundary(:,[1,4,3,2]);
    o(:,:,2) = boundary(:,[1,2,6,5]);
    o(:,:,3) = boundary(:,[2,3,7,6]);
    o(:,:,4) = boundary(:,[3,4,8,7]);
    o(:,:,5) = boundary(:,[4,1,5,8]);
    o(:,:,6) = boundary(:,[5,6,7,8]);
    for i = 1:6
        p=patch('Vertices',nodes,'Faces',s(:,:,i));
        set(p,'facecolor','yellow','edgecolor','black');
    end
    for i = 1:6
        p=patch('Vertices',nodes,'Faces',o(:,:,i));
        set(p,'facecolor','red','edgecolor','black');
    end
end

view(50,50)

end

