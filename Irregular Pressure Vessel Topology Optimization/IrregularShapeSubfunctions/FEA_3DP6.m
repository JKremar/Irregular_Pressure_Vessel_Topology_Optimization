function [U,K,F] = FEA_3DP6(struc,elements,map,KE,Po,noF,fixeddofs,oldstruc,oldK,oldF)
%for irregular shapes

%Computes Finite element analysis for the structure where 1 means there is
%material and 0 corresponds to void
%Inputs:    -struc:material distribution representation (1=material & 0=void)
%           -elements:mapping of which nodes belong to each element and
%           their relative positionings
%           -KE:elemental k matrix
%           -Po:magnitude of pressure
%           -e:used in computation of diriac delta function to determine
%                   pressure loading for the given LSF
%           -fixeddofs:[degrees of freedom that are fixed]
%           -oldstruc:The previous structure that the K matrix was
%                   calculated for so that the elements that don't change don't
%                   need to be recomputed in the global K-matrix
%           -oldK:Previous global K matrix to serve as starting point for
%                   the this iteration
%Outputs:   -U:dispacement vector result from FEA
%           -K:current global K matrix to be used as starting point for the
%                   next iteration



%Initialize F,K and U Matrices--------------------------------
numnodes=max(max(elements));
numelements=size(elements,1);
U=zeros(3*numnodes,1);
dof=3*repelem(elements,1,3)-repmat([2,1,0],1,8);
fe=repmat(Po',8,1).*[1;1;1;-1;1;1;-1;-1;1;1;-1;1;1;1;-1;-1;1;-1;-1;-1;-1;1;-1;-1];
%-------------------------------------------------------------

%Compute pressure force---------------------------------------

%-------------------------------------------------------------

if(nargin<7 || isempty(oldK)|| isempty(oldF))
    %compute full K matrix
    F=zeros(3*numnodes,1);
    K=sparse(3*numnodes,3*numnodes);
    %eKE=permute(max(struc(map),0.0001),[3,2,1]).*KE;
    for(e=1:numelements)
        %K(dof(e,:),dof(e,:))=K(dof(e,:),dof(e,:))+eKE(:,:,e);
        K(dof(e,:),dof(e,:))=K(dof(e,:),dof(e,:))+max(struc(map(e)),0.0001)*KE;
        if(struc(map(e))==1)
            F(dof(e,:))=F(dof(e,:))+fe;
        end
%         if(~mod(e,500))
%             fprintf('Assembled %d elements of %d \n',e,numelements);
%         end
    end
else
    %only modify K and F where needed
    K=oldK;  F=oldF;
    ele=find(struc(map)-oldstruc(map));    %elements that changed
    for(i=1:numel(ele))     %0==no change, 1==added material, -1==removed material
        Ke_old=max(oldstruc(map(ele(i))),0.0001)*KE;
        Ke=max(struc(map(ele(i))),0.0001)*KE;
        K(dof(ele(i),:),dof(ele(i),:))=K(dof(ele(i),:),dof(ele(i),:))-Ke_old+Ke;
        Fe_old=oldstruc(map(ele(i)))*fe;
        Fe=struc(map(ele(i)))*fe;
        F(dof(ele(i),:))=F(dof(ele(i),:))-Fe_old+Fe;
    end
    
end

%Solve System of Equations
F(noF)=0;
freedofs=setdiff(1:3*numnodes,fixeddofs);
U(freedofs,:)=K(freedofs,freedofs)\F(freedofs,:);


end
