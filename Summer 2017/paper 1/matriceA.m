function [ A] = matriceA(kt,C,G,D,a,h)
W=sqrt(9.81/(C*h));
N=9;%change to 9
i=sqrt(-1);
k=kt;
  syms x
for p=1:N+1
    for q=1:N+1
        if p==1 & q==1
            a1(p,q)=int(conj(cosh(k(q)*(x+h))/cosh(k(q)*h))*(cosh(k(p)*(x+h))/cosh(k(p)*h)),[-h -a]);
        elseif p==1 & q~=1
            a1(p,q)=int(conj(cos(k(q)*(x+h))/cos(k(q)*h))*(cosh(k(p)*(x+h))/cosh(k(p)*h)),[-h -a]);
        elseif p~=1 & q==1
            a1(p,q)=int(conj(cosh(k(q)*(x+h))/cosh(k(q)*h))*(cos(k(p)*(x+h))/cos(k(p)*h)),[-h -a]);
        else
            a1(p,q)=int(conj(cos(k(q)*(x+h))/cos(k(q)*h))*(cos(k(p)*(x+h))/cos(k(p)*h)),[-h -a]);
        end
    end
end
syms x
%for -a<= y <=0
C1=i*(cosh(k(1)*(x+h))/cosh(k(1)*h))*(k(1)+2*(k(1)*G+(i*W^2/(D*k(1)^2)))) - (W^2/D)*((x/a)*((2*cosh(k(1)*(h-a))/(k(1)^2*cosh(k(1)*h)))-(2/k(1)^2)) - 2/k(1)^2);
syms y x
c1=@(y) -1*(cos(y*(x+h))/cos(y*h));
c2=@(y)   y-2*i*k(1)*G -2*W^2/(D*y^2);
c3=@(y)   -(x/a)*((2*cos(y*(h-a))/(y^2*cos(y*h)))+(2/y^2)) + 2/y^2 ;
Cn=@(y) c1(y)*c2(y) - (W^2/D)*c3(y);
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

