function [ R ] = scatter( i,T,X,G,C,h,theta,N)

syms k0i
        sol=solve(1==k0i*C*tanh(k0i),k0i);
        k0=abs(sol);
syms x
limcap=1.57;
%it can be observed roots are always between npi - atan(c)  and npi
 for p=1:N
F=@(K) C*K+cot(K);
X=[double((p*pi-(pi-limcap))) double((p*pi-(pi-limcap)+1.57))];
z(p) = fzero(F,X);
 end
 %now z is the array of k1*h k2*h .....
    for v=1:N+1
        if v==1
            kt(v)=k0;
        else
            kt(v)=z(v-1);
        end
    end
 kt=double(kt*(1/h));

k=kt;




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
ri=linsolve(A,B);
R=ri(1);
end
