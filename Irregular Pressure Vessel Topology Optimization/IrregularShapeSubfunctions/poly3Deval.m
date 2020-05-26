function [eval] = poly3Deval(coef,points)
%

if(iscell(coef))
    [Req,Ceq]=size(coef);
else
    Req=1;
    Ceq=1;
end
if(nargin==1)
    P=input('at what location would you like to evaluate? :');
else
    P=points;
end

eval=zeros([Req,Ceq,size(points,2)]);         %Initializes function value matrix
for(r=1:Req)
    for(c=1:Ceq)
        for(p=1:size(points,2))
            xp=0;
            for(i=size(coef{r,c},1):-1:1)
                yp=0;
                for(j=size(coef{r,c},2):-1:1)
                    zp=0;
                    for(k=size(coef{r,c},3):-1:1)
                        eval(r,c,p)=eval(r,c,p)+coef{r,c}(i,j,k)*P(1,p)^xp*P(2,p)^yp*P(3,p)^zp;
                        zp=zp+1;
                    end
                    yp=yp+1;
                end
                xp=xp+1;
            end
        end
    end
end

end

