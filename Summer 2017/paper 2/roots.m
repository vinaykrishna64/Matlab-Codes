function [ kt] = roots(C,h,N)

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



end

