function [ Xf ] = chanleeextended(C,G,h,a,H,mn,tn  )

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
%calculation of constants
i=sqrt(-1);
W=sqrt(9.81/(C*h));
rho=1025;
g=9.81;
T=tn*rho*g*h^2;
m=mn*rho*h;
phi_0=-i*g*H/(2*W);
P=i*sqrt(m*W^2/T);

%N=9;%change to 9

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
temp1= i*k(1)*phi_0*(1+2*G)*(cosh(k(1)*(x+h))/cosh(k(1)*h));
temp2= (exp(P*a)-(cosh(k(1)*(h-a))/cosh(k(1)*h)))*(sinh(P*x)/sinh(-P*a));
temp3= (cosh(k(1)*(x+h))/cosh(k(1)*h))-exp(-P*x);
A_0= -rho*g*h/(T*k(1)^2);
C1=temp1 - i*W*A_0*(temp2+temp3);


syms y x
temp11=@(y) -1*phi_0*(y-2*i*G*k(1))*(cos(y*(x+h))/cos(y*h));
temp21= @(y)(exp(p*a)-(cos(y*(h-a))/cos(y*h)))*(sinh(P*x)/sinh(-P*a));
temp31= @(y)(cos(y*(x+h))/cos(y*h))-exp(-P*x);
A_n=@(y) 2*i*rho*W/(T*y^2);
Cn=@(y) temp11(y)-i*W*A_n(y)*(temp21(y)+temp31(y));
Cnst=@(y) conj(temp11(y)-i*W*A_n(y)*(temp21(y)+temp31(y)));
%the above are the coefficients of f0 and fn's in the barrier region
syms x
for p=1:N+1
    for q=1:N+1
        if p==1 & q==1
            a2(p,q)=int(conj(C1)*C1,[-a 0]);
        elseif p==1 & q~=1
            a2(p,q)=int(Cnst(k(q))*C1,[-a 0]);
        elseif p~=1 & q==1
            a2(p,q)=int(conj(C1)*Cn(k(p)),[-a 0]);
        else
            a2(p,q)=int(Cnst(k(q))*Cn(k(p)),[-a 0]);
        end
    end
end

A=a1+a2;
A=double(A); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x
for p=1:N+1
    
    if p==1
        B(p)=-1*int(conj(-i*phi_0*(k(1)/cosh(k(1)*h))*cosh(k(1)*(x+h)))*C1,[-a 0]);
    else
        B(p)=-1*int(conj(-i*phi_0*(k(1)/cosh(k(1)*h))*cosh(k(1)*(x+h)))*Cn(k(p)),[-a 0]);
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

