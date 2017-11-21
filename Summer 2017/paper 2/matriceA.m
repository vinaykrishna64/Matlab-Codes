function [ A ] = matriceA(i,T,X,N,k,G,theta)
%pn=sqrt(kn^2+v^2)
%e=k0*cos(theta) v=k0 sintheta
e=k(1)*cosd(theta);
v=k(1)*sind(theta);
ci=sqrt(-1);
for p=1:N+1
    for q=1:N+1
       A(p,q)=0;
    end
end
l=1;
const1=ci*(e+2*k(1)*G);
const2=2*ci*G*k(1);
syms y
while l<i+1
for p=1:N+1
    p1p=sqrt(k(p)^2+v^2);
    for q=1:N+1
        if T(l)=='G'
            if p==1 & q==1
             A(p,q)=int(conj(cosh(k(1)*y))*cosh(k(1)*y),X{l})+A(p,q);
            elseif p==1 & q~=1
                A(p,q)=int(conj(cos(k(q)*y))*cosh(k(1)*y),X{l})+A(p,q);
            elseif p~=1 & q==1
                A(p,q)=int(conj(cosh(k(1)*y))*cos(k(p)*y),X{l})+A(p,q);
            else
                A(p,q)=int(conj(cos(k(q)*y))*cos(k(p)*y),X{l})+A(p,q);
            end
        elseif T(l)=='B'
            p1=sqrt(k(q)^2+v^2);
            
             if p==1 & q==1
             A(p,q)=int(conj(const1*cosh(k(1)*y))*const1*cosh(k(1)*y),X{l})+A(p,q);
            elseif p==1 & q~=1
                A(p,q)=int(conj(-(p1-const2)*cos(y*k(q)))*const1*cosh(k(1)*y),X{l})+A(p,q);
            elseif p~=1 & q==1
                A(p,q)=int(conj(const1*cosh(k(1)*y))*-(p1p-const2)*cos(y*k(p)),X{l})+A(p,q);
             else
                A(p,q)=int(conj(-(p1-const2)*cos(y*k(q)))*-(p1p-const2)*cos(y*k(p)),X{l})+A(p,q);
            end
    end
    end
end
l=l+1;

end
end
