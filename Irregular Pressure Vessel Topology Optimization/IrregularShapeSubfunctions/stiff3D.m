function [Ke,B,C] = stiff3D(E,v,lx,ly,lz)
%Calculates the elemental stiffness matrices

%Inputs:    -E:modulus of elasticity
%           -v:Poison's ratio
%Outputs:   -Ke:elemental stiffness matrix
%           -Ktr:element matrix for trace tensor
%           -lamda:Lame Constant 

C=(E/((1+v)*(1-2*v)))*[[(1-v)*eye(3)+v*~eye(3)],zeros(3);zeros(3),((1-2*v)/2)*eye(3)];



if(nargin==2)
    %ke=[(3-v)/6 , (1+v)/8 , (-3-v)/12 , (3*v-1)/8 , (v-3)/12 , (-1-v)/8 , v/6 , (1-3*v)/8];
    kp=[-(3*v-2)/9,1/24,-1/18,-(4*v-1)/24,(4*v-1)/24,1/36,1/48,-1/24,(6*v-5)/72,-(4*v-1)/48,-1/48,(4*v-1)/48,(3*v-1)/36,(3*v-2)/36];
    k1=kp([1,2,2,3,5,5;2,1,2,4,6,7;2,2,1,4,7,6;3,4,4,1,8,8;5,6,7,8,1,2;5,7,6,8,2,1]);
    k2=kp([9,8,12,6,4,7;8,9,12,5,3,5;10,10,13,7,4,6;6,5,11,9,2,10;4,3,5,2,9,12;11,4,6,12,10,13]);
    k3=kp([6,7,4,9,12,8;7,6,4,10,13,10;5,5,3,8,12,9;9,10,2,6,11,5;12,13,10,11,6,4;2,12,9,4,5,3]);
    k4=kp([14,11,11,13,10,10;11,14,11,12,9,8;11,11,14,12,8,9;13,12,12,14,7,7;10,9,8,7,14,11;10,8,9,7,11,14]);
    k5=kp([1,2,8,3,5,4;2,1,8,4,6,11;8,8,1,5,11,6;3,4,5,1,8,2;5,6,11,8,1,8;4,11,6,2,8,1]);
    k6=kp([14,11,7,13,10,12;11,14,7,12,9,2;7,7,14,10,2,9;13,12,10,14,7,11;10,9,2,7,14,7;12,2,9,11,7,14]);

    Ke=(E/((1+v)*(1-2*v)))*[k1,k2,k3,k4;k2',k5,k6,k3';k3',k6,k5',k2';k4,k3,k2,k1'];

    dN_cent=0.25*[0,-1,-1,-1;0,1,-1,-1;0,1,1,-1;0,-1,1,-1;0,-1,-1,1;0,1,-1,1;0,1,1,1;0,-1,1,1]';
    order=[1,0,0;0,2,0;0,0,3;0,3,2;3,0,1;2,1,0];
    B=dN_cent([order+1,order+5,order+9,order+13,order+17,order+21,order+25,order+29]);
else
    if(nargin==3)
        if(numel(lx)==3)
            ly=lx(2);
            lz=lx(3);
            lx=lx(1);
        else
            ly=lx(1);
            lz=lx(1);
            lx=lx(1);
        end
    end
    num_nodes=8;
    J=[lx/2,ly/2,lz/2].*eye(3);
    dN=cell(8,3);
    n=[-1/2,1/2;1/2,1/2];   dn=[-1/2;1/2];  %1-D shape function coefficients and derivative coefficients
    for(i=1:num_nodes)
        xy=n(floor(mod(i-1,4)/2)+1,:)'*n(floor(mod(i,4)/2)+1,:);     %[y]'*[x]
        dxy=n(floor(mod(i-1,4)/2)+1,:)'*dn(floor(mod(i,4)/2)+1,:);   %[y]'*[dx]
        xdy=dn(floor(mod(i-1,4)/2)+1,:)'*n(floor(mod(i,4)/2)+1,:);   %[dy]'*[x]
        for(c=1:size(xy,2))
            if(c<=size(dxy,2))      %because partial wrt x will have 1 less column
                dN{i,1}=[dN{i,1},permute(dxy(:,c)*n(floor((i-1)/4)+1,:),[1,3,2])];    %[dxy(:,c)]*[z]    
            end
            dN{i,2}=[dN{i,2},permute(xdy(:,c)*n(floor((i-1)/4)+1,:),[1,3,2])];   %[xdy(:,c)]*[z]
            dN{i,3}=[dN{i,3},permute(xy(:,c)*dn(floor((i-1)/4)+1,:),[1,3,2])];   %[xy(:,c)]*[dz]
        end
    end

    P_1D=[-1/3^0.5,1/3^0.5];    W=ones(1,8);
    GPts=P_1D([1,1,1,1,2,2,2,2;1,1,2,2,1,1,2,2;1,2,1,2,1,2,1,2]);
    dNg=poly3Deval(dN,GPts);
    [num_nodes,D,G]=size(dNg);
    Ke=zeros(num_nodes*D);
    for(g=1:G)
        delNg=J^-1*dNg(:,:,g)';
        B=zeros(2*D,num_nodes*D);
        for(n=1:num_nodes)
            for(d=1:D)
                a=[1:d-1,d+1:3];
                B(d,n*D-D+d)=delNg(d,n);
                B(D+a(1),n*D-D+d)=delNg(a(2),n);
                B(D+a(2),n*D-D+d)=delNg(a(1),n);
            end
        end
        Ke=Ke+B'*C*B*det(J)*W(g);
    end
    dN0=poly3Deval(dN,[0;0;0]);
    delN0=J^-1*dN0';
    B=zeros(2*D,num_nodes*D);
    for(n=1:num_nodes)
        for(d=1:D)
            a=[1:d-1,d+1:3];
            B(d,n*D-D+d)=delN0(d,n);
            B(D+a(1),n*D-D+d)=delN0(a(2),n);
            B(D+a(2),n*D-D+d)=delN0(a(1),n);
        end
    end
end



end

