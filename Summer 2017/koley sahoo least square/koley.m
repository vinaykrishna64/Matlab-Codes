function [R0] = koley(h,w,G,theta,Rm,RT,a,b)

%input required
%a,b for length of membrane
%h,w,G,theta,Rm,RT
g=9.801; %acceleration due to gravity
H=h/10;%approximation in the paper
pi=3.1416;
rho=1025;
i=sqrt(-1);

%%%% Assuming N=6
N=6;%total number of evanascent wave numbers
syms k0 kn
k0=solve(w^2==k0*g*tanh(k0*h),k0);
k(1)=abs(k0);
C=w^2*h/g;
for p=1:N
    X=[(p-1)*pi+(pi-1.4) p*pi];
    f=@(kn) tan(kn)+C/kn;
    k(p+1)=fzero(f,X)/h;
end
k=double(k);
N=N+1;

%all coefficients
m=rho*h*Rm;%mass
T=RT*rho*g*h^2;%tension
beta0=k(1)*sind(theta);
alpha0=k(1)*cosd(theta);
alphan=@(y) sqrt(beta0^2+y^2);
a0=-rho*g*H/(T*(k(1)^2+(m*w^2/T)-beta0^2));
an=@(y) 2*i*rho*w/(T*(y^2-(m*w^2/T)-beta0^2));
g0=alpha0*g*H/(2*w);
b0=i*w*a0-(alpha0*g0*H/(2*w));
bn=@(y) i*w*an(y)+alphan(y);
c0=k(1)*g*H*G/w;
cn=2*i*k(1)*G;
delta=i*sqrt((m*w^2/T)-beta0^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms x
%%%%%%%%%%%%%      MATRIX SYSTEMS
%A*[A B A0 A1 .....]'=B
%% A
for hor=1:N+2
    for ver=1:N+2
         %the two equations from barrier displacement condition 
       if hor<3
              A(hor,ver)=0;
            
        %the system from potential equations
       
            elseif ver==1 & hor==3
                A(hor,ver)=i*w*int(exp(delta*x)*(cosh(k(1)*(h-x))/cosh(k(1)*h)),[a b]);
            elseif ver==2 & hor==3
                A(hor,ver)=i*w*int(exp(-delta*x)*(cosh(k(1)*(h-x))/cosh(k(1)*h)),[a b]);
            elseif ver==3 & hor==3
                A(hor,ver)=b0*int((cosh(k(1)*(h-x))/cosh(k(1)*h))^2,[a b]);
            elseif ver>3 & hor==3
                A(hor,ver)=bn(k(ver-2))*int((cos(k(ver-2)*(h-x))/cos(k(ver-2)*h))*(cosh(k(1)*(h-x))/cosh(k(1)*h)),[a b]);    
            elseif ver==1 & hor>3
                A(hor,ver)=i*w*int(exp(delta*x)*(cos(k(hor-2)*(h-x))/cos(k(hor-2)*h)),[a b]);
            elseif ver==2 & hor>3
                A(hor,ver)=i*w*int(exp(-delta*x)*(cos(k(hor-2)*(h-x))/cos(k(hor-2)*h)),[a b]);
            elseif ver==3 & hor>3
                A(hor,ver)=b0*int((cos(k(hor-2)*(h-x))/cos(k(hor-2)*h))*(cosh(k(1)*(h-x))/cosh(k(1)*h)),[a b]);
            elseif ver>3 & hor>3
                A(hor,ver)=bn(k(ver-2))*int((cos(k(ver-2)*(h-x))/cos(k(ver-2)*h))*(cos(k(hor-2)*(h-x))/cos(k(hor-2)*h)),[a b]);
            end
            
        end
        
    end


for hor=1:N+2
    for ver=1:N+2
     %the two equations from barrier displacement condition 
    if hor<3
            if ver==1 & hor==1
                A2(hor,ver)=exp(delta*a);
            elseif ver==2 & hor==1
                A2(hor,ver)=exp(-delta*a);
            elseif ver==3 & hor==1
                A2(hor,ver)=a0*cosh(k(1)*(h-a))/cosh(k(1)*h);
            elseif ver>3 & hor==1
                A2(hor,ver)=an(k(ver-2))*cos(k(ver-2)*(h-a))/cos(k(ver-2)*h);    
            elseif ver==1 & hor==2
                A2(hor,ver)=exp(delta*b);
            elseif ver==2 & hor==2
                A2(hor,ver)=exp(-delta*b);
            elseif ver==3 & hor==2
                A2(hor,ver)=a0*cosh(k(1)*(h-b))/cosh(k(1)*h);
            elseif ver>3 & hor==2
              
               A2(hor,ver)=an(k(ver-2))*cos(k(ver-2)*(h-b))/cos(k(ver-2)*h);
            end
      else
    
            if hor==ver & hor>2
                if hor==3
                    A2(hor,ver)=c0*int((cosh(k(1)*(h-x))/cosh(k(1)*h))^2,[0 h]);
                else
                    A2(hor,ver)=cn*int((cos(k(hor-2)*(h-x))/cos(k(hor-2)*h))^2,[0 h]);
                end
            else
                A2(hor,ver)=0;
            end
        
        end
    end
end



%% B
for hor=1:N+2
    syms x
    if hor<3
        B(hor)=0;
    elseif hor==3
    B(hor)=int(g0*(cosh(k(1)*(h-x))/cosh(k(1)*h))^2,[a,b]);
    else
    B(hor)=int(g0*(cosh(k(1)*(h-x))/cosh(k(1)*h))*(cos(k(hor-2)*(h-x))/cos(k(hor-2)*h)),[a,b]);        
    end
end
if iscolumn(B)==0
    B=B';
end
R=linsolve(A2-A,B);
R0=double(abs(R(3)));
end