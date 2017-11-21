function [ B] = matriceB( kt,C,G,D,a,h )
W=sqrt(9.81/(C*h));
N=9;%change to 9
i=(-1)^(0.5);
k=kt;
syms x
%for -a<= y <=0
C1=i*(cosh(k(1)*(x+h))/cosh(k(1)*h))*(k(1)+2*(k(1)*G+(i*W^2/(D*k(1)^2)))) - (W^2/D)*((x/a)*((2*cosh(k(1)*(h-a))/(k(1)^2*cosh(k(1)*h)))-(2/k(1)^2)) - 2/k(1)^2);
syms y x
Cn=@(y) -1*(cos(y*(x+h))/cos(y*h))*(y-2*i*(k(1)*G-(i*W^2/(D*y^2)))) - (W^2/D)*(-(x/a)*((2*cos(y*(h-a))/(y^2*cos(y*h)))+(2/y^2)) + 2/y^2);
%the above are the coefficients of f0 and fn's in the barrier region
syms x
for p=1:N+1
    if p==1
        B(p)=-1*int(conj(-i*(k(1)/cosh(k(1)*h))*cosh(k(1)*(x+h)))*C1,[-a 0]);
    else
        B(p)=-1*int(conj(-i*(k(1)/cosh(k(1)*h))*cosh(k(1)*(x+h)))*Cn(k(p)),[-a 0]);
end

B=double(B);
if iscolumn(B)==0
    B=B';
end

end

