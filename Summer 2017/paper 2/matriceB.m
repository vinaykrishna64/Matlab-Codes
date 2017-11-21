function [ B ] = matriceB( i,T,X,N,k,G,theta )
%pn=sqrt(kn^2+v^2)
%e=k0*cos(theta) v=k0 sintheta
e=k(1)*cosd(theta);
v=k(1)*sind(theta);
ci=sqrt(-1);
for p=1:N+1
   B(p)=0; 
end
l=1;
syms y
while l<i+1
for p=1:N+1
    if T(l)=='G'
            if p==1 
                B(p)=0+B(p);
            else
                 B(p)=0+B(p);
            end
    elseif T(l)=='B'
             if p==1 
                 B(p)=int(conj(-ci*e*cosh(k(1)*y))*ci*(e+2*k(1)*G)*cosh(k(1)*y),X{l})+B(p);
            else
                B(p)=int(conj(-ci*e*cosh(k(1)*y))*-(sqrt(k(p)^2+v^2)-(2*ci*G*k(1)))*cos(y*k(p)),X{l})+B(p);
             end
    end

end
l=l+1;

end
if iscolumn(B)==0
    B=B';
end
end

