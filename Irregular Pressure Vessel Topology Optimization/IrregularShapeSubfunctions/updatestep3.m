function [lsf] = updatestep3(lsf,shapeSens,stepLength,bearing,Le)
%updates the structure and the level-set function

%smooth sensitivities
%C=reshape([0,0,0;0,1,0;0,0,0;0,1,0;1,2,1;0,1,0;0,0,0;0,1,0;0,0,0],3,3,3)/8;
C=reshape([0,1,0;1,2,1;0,1,0;1,2,1;2,3,2;1,2,1;0,1,0;1,2,1;0,1,0],3,3,3)/27;
shapeSens=convn(padarray(shapeSens,[1,1,1],'replicate'),C,'valid');

%Insure load bearing pixels remain solid
shapeSens(bearing)=0;

v=-shapeSens;
%add zeros to boarder of v 
vFull=zeros(size(v)+2); vFull(2:end-1,2:end-1,2:end-1)=v;
lsf=padarray(lsf,[1,1,1],'replicate');

%determine timestep (based on CFL condition)
dt=Le*0.1/max(abs(v(:)));

for(i=1:(10*stepLength))
    dpx=circshift(lsf,[-1,0,0])-lsf;  %Find derivatives on the grid
    dmx=lsf-circshift(lsf,[1,0,0]);
    dpy=circshift(lsf,[0,-1,0])-lsf;
    dmy=lsf-circshift(lsf,[0,1,0]);
    dpz=circshift(lsf,[0,0,-1])-lsf;
    dmz=lsf-circshift(lsf,[0,0,1]);
    %Update LSF
    lsf=lsf-dt*min(vFull,0).*sqrt(min(dmx,0).^2+max(dpx,0).^2+min(dmy,0).^2+max(dpy,0).^2+min(dmz,0).^2+max(dpz,0).^2) ...
        -dt*max(vFull,0).*sqrt(max(dmx,0).^2+min(dpx,0).^2+max(dmy,0).^2+min(dpy,0).^2+max(dmz,0).^2+min(dpz,0).^2);
end

lsf=lsf(2:end-1,2:end-1,2:end-1);

end

