function [ Xf ] = chanlee(C,G,D,a,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=9;%nine
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
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W=sqrt(9.81/(C*h));
%N=9;%change to 9
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
coeff=(i/cosh(k(1)*h)) *(k(1) + 2*G*k(1) + (2i*W^2 / (D* k(1)^2)));
coeff2=((2/k(1)^2)*(cosh(k(1)*(h-a))/cosh(k(1)*h))-(2/k(1)^2))/a;
C1=coeff*cosh(k(1)*(h+x)) -(W^2/D)*(coeff2*x-(2/k(1)^2));
syms y x
c1=@(y) -1*(cos(y*(x+h))/cos(y*h));
c2=@(y)   y-2*i*k(1)*G -(2*W^2/(D*y^2));
c3=@(y)   -(x/a)*((2*cos(y*(h-a))/(y^2*cos(y*h)))-(2/y^2)) + 2/y^2 ;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x
for p=1:N+1
    
    if p==1
        B(p)=-1*int(conj(-i*(k(1)/cosh(k(1)*h))*cosh(k(1)*(x+h)))*C1,[-a 0]);
    else
        B(p)=-1*int(conj(-i*(k(1)/cosh(k(1)*h))*cosh(k(1)*(x+h)))*Cn(k(p)),[-a 0]);
end
end
B=double(B);
if iscolumn(B)==0
    B=B';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Xi=linsolve(A,B);
Xf=abs(Xi(1));

end

