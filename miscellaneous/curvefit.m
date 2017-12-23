function [ Y_new,coeff ] = curvefit( x,y,N)
%curvefit( X,Y,N)
%takes X and Y values and fits nth order curve
%also gives new values and coefficient of the polynomial 
if isrow(x)==0
    x=x';
end
if isrow(y)==0
    y=y';
end
x1=ones(length(x),N+1);
coeff=zeros(N+1,1);
for i=2:N+1
    x1(:,i)=(x.^(i-1))';
end
coeff=inv(x1'*x1)*x1'*y';
Y_new=x1*coeff;

end

