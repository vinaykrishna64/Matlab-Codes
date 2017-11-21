function [ B] = newB( kt,C,G,D,a,h,angle )
W=sqrt(9.81/(C*h));
N=9;%change to 9
i=sqrt(-1);
k=kt;
E=k(1)*cosd(angle);
V=k(1)*sind(angle);
syms x
%for -a<= y <=0
C1=i*cosh(k(1)*(x+h))*(E+2*(k(1)*G+(i*W^2/(D*k(1)^2)))) - (W^2/D)*((x/a)*((2/k(1)^2)*cosh(k(1)*(h-a))-(2/k(1)^2)) -(2/k(1)^2));
syms y x
Cn=@(y) -1*cos(y*(x+h))*(sqrt(y^2+V^2)-2*i*(k(1)*G-(i*W^2/(D*y^2)))) - (W^2/D)*((x/a)*(-(2/y^2)*cos(y*(h-a))+(2/y^2)) +(2/y^2));
%the above are the coefficients of f0 and fn's in the barrier region
syms x
for p=1:N+1
    if p==1
        B(p)=-1*int(conj(-i*E*cosh(k(1)*(x+h)))*C1,[-a 0]);
    else
        B(p)=-1*int(conj(-i*E*cosh(k(1)*(x+h)))*Cn(k(p)),[-a 0]);
end

B=double(B);
if iscolumn(B)==0
    B=B';
end

end

