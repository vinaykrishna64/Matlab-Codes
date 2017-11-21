function [ A] = newA(kt,C,G,D,a,h,angle)
W=sqrt(9.81/(C*h));
N=9;%change to 9
i=sqrt(-1);
k=kt;
E=k(1)*cosd(angle);
V=k(1)*sind(angle);
  syms x
for p=1:N+1
    for q=1:N+1
        if p==1 & q==1
            a1(p,q)=int(conj(cosh(k(q)*(x+h)))*(cosh(k(p)*(x+h))),[-h -a]);
        elseif p==1 & q~=1
            a1(p,q)=int(conj(cos(k(q)*(x+h)))*(cosh(k(p)*(x+h))),[-h -a]);
        elseif p~=1 & q==1
            a1(p,q)=int(conj(cosh(k(q)*(x+h)))*(cos(k(p)*(x+h))),[-h -a]);
        else
            a1(p,q)=int(conj(cos(k(q)*(x+h)))*(cos(k(p)*(x+h))),[-h -a]);
        end
    end
end
syms x
%for -a<= y <=0
C1=i*cosh(k(1)*(x+h))*(E+2*(k(1)*G+(i*W^2/(D*k(1)^2)))) - (W^2/D)*((x/a)*((2/k(1)^2)*cosh(k(1)*(h-a))-(2/k(1)^2)) -(2/k(1)^2));
syms y x
Cn=@(y) -1*cos(y*(x+h))*(sqrt(y^2+V^2)-2*i*(k(1)*G-(i*W^2/(D*y^2)))) - (W^2/D)*((x/a)*(-(2/y^2)*cos(y*(h-a))+(2/y^2)) +(2/y^2));
%the above are the coefficients of f0 and fn's in the barrier region
syms x
for p=1:N+1
    for q=1:N+1
        if p==1 & q==1
            a2(p,q)=int(conj(C1)*C1,[-a 0]);
        elseif p==1 & q~=1
            a2(p,q)=int(conj(Cn(k(q)))*C1,[-a 0]);
        elseif p~=1 & q==1
            a2(p,q)=int(conj(C1)*Cn(k(p)),[-a 0]);
        else
            a2(p,q)=int(conj(Cn(k(q)))*Cn(k(p)),[-a 0]);
        end
    end
end


A=a1+a2;
A=double(A);   
end

