function [z] = poly2Deval(coef,x,y)
%finds the function values of the polynomial with coefficients 'coef' evaluated
%at x and y

%inputs:    coef:coefficient matrix of the polynomial where
%               -bottom right value is the constant
%               -moving left increases the power of x
%               -moving up increases the power of y
%           x:x-values to evaluate the function at
%           y:y-values to evaluate the function at

%outputs:   z: the values of the function evaluated at the x values

if(nargin==1)
    P=input('at what location would you like to evaluate? :');
elseif(nargin==2)            %If the user imputted the point as a vector
    P(:,1,1)=x(:,1);
    P(:,1,2)=x(:,2);
else
    for(i=1:numel(x))
        for(j=1:numel(y))
            P(i,j,:)=[x(i),y(j)];
        end
    end
end
    

z=zeros(size(P,1),size(P,2));         %Initializes function value matrix
for(i=1:size(P,1))
    for(j=1:size(P,2))
        b=0;
        for(r=size(coef,1):-1:1)
            a=0;
            for(c=size(coef,2):-1:1)
                z(i,j)=z(i,j)+coef(r,c)*P(i,j,1)^a*P(i,j,2)^b;  %Adds terms to z
                a=a+1;                                  %Increases power of y
            end
            b=b+1;                                      %Increases power of x
        end
    end
end


end

