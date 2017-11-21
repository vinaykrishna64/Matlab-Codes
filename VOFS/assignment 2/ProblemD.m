clear all
clc
tic
rho=7850;
E=209*10^9;
L=100; 
B=2;
h=1;
A=B*h;
I=(1/12)*B*h^3;
k=1000;%ktheta
xv=linspace(0,100,100);
syms b
mat=[[-1 0 1 0];
    [0 -1 0 1];
    [cos(b) sin(b) cosh(b) sinh(b)];
    [(-b*cos(b)+k*sin(b)) (-b*sin(b)-k*cos(b)) (b*cosh(b)-k*sinh(b)) (b*sinh(b)-k*cosh(b))]];
f=det(mat);
fd=diff(f);
f=matlabFunction(f);
fd=matlabFunction(fd);

for j=1:10
t=2+(j-1)*pi;
for i=1:500
    t=t-f(t)/fd(t);
end
bl(j)=t;
end
syms g2 
mat1=mat*[1;g2;1;g2];

fn=matlabFunction(mat1(3));
for j=1:10
    s=solve([fn(bl(j),g2)==0],[g2]);
    g2v(j)=s;
end

for j=1:10
    for i=1:length(xv)
        phi(j,i)=cos(bl(j)*xv(i)/L)+g2v(j)*sin(bl(j)*xv(i)/L)+cosh(bl(j)*xv(i)/L)+g2v(j)*sinh(bl(j)*xv(i)/L);
    end
end

plot(xv,phi)
for j=1:10
    for k=1:10
        phi2(j,k)=0;
        for i=1:length(xv)
            phi2(j,k)=phi2(j,k)+phi(j,i).*phi(k,i)*(xv(2)-xv(1));
        end
    end
end
figure
mesh(phi2)


%%%% mode summation %%%%%%%%%%%%%%
for i=1:10
    for j=1:10
 if i==j
     mass(i,j)=rho*A*phi2(i,i);
 else
     mass(i,j)=0;
 end
    end
end
figure
mesh(mass)

for i=1:10
    for j=1:10
 if i==j
     stiff(i,j)=E*I*(bl(j)/L)^4*phi2(i,i);
 else
     stiff(i,j)=0;
 end
    end
end
figure
mesh(stiff)

for i=1:10
    for j=1:10
        if i==j
            damp(i,j)=2*0.05*sqrt(mass(i,i)*stiff(i,i));
        else
            damp(i,j)=0;
        end
    end
end
figure
mesh(damp)
M1=mass;
M2=damp;
M3=stiff;
zeta=M2(1,1)/(2*sqrt(M1(1,1)*M3(1,1)));
%%%FK force
w=2;
k=w^2/9.81;
T=.5;
a=1/(7*k);
syms x t
fk=a*cos(k*x-2*t)*exp(k*T)*9.81*1.025*B;
figure
fsurf(a*cos(k*x-2*t)*exp(k*T)*9.81*1.025*B)
alp1=atan2(-w*M2(1,1),(M3(1,1)-M1(1,1)*w^2));
fk=a*cos(k*x-2*(t-alp1))*exp(k*T)*9.81*1.025*B;
syms t
PI=real(int(fk*(cos(bl(1)*x/L)+g2v(1)*sin(bl(1)*x/L)+cosh(bl(1)*x/L)+g2v(1)*sinh(bl(1)*x/L)),x,0,L)/sqrt((M3(1,1)-M1(1,1)*w^2)^2+(M2(1,1)*w)^2));
CF=real(exp(-zeta*sqrt(M3(1,1)/M1(1,1))*t)*cos(sqrt(M3(1,1)/M1(1,1))*sqrt(1-zeta^2)*t));
total=matlabFunction(CF+PI);
CF=matlabFunction(CF);
figure
fplot(CF,[0,100])
PI=matlabFunction(PI);
figure
fplot(PI,[0,100])
figure
fplot(total,[0,100])

for j=1:10
    alp=atan2(-w*M2(j,j),(M3(j,j)-M1(j,j)*w^2));
    fk=a*cos(k*x-2*(t-alp))*exp(k*T)*9.81*1.025*B;
    qm(j)=int((cos(bl(j)*x/L)+g2v(j)*sin(bl(j)*x/L)+cosh(bl(j)*x/L)+g2v(j)*sinh(bl(j)*x/L))*fk,x,0,L)/(stiff(j,j)-4*mass(j,j));
end
syms t
qm(4)=cos(sqrt(stiff(1,1)/mass(1,1))*t)+qm(4);
for j=1:10
    phim(j)=cos(bl(j)*x/L)+g2v(j)*sin(bl(j)*x/L)+cosh(bl(j)*x/L)+g2v(j)*sinh(bl(j)*x/L);
end
for j=1:10
    phid(j)=diff(cos(bl(j)*x/L)+g2v(j)*sin(bl(j)*x/L)+cosh(bl(j)*x/L)+g2v(j)*sinh(bl(j)*x/L));
end
for j=1:10
    phi3d(j)=diff(cos(bl(j)*x/L)+g2v(j)*sin(bl(j)*x/L)+cosh(bl(j)*x/L)+g2v(j)*sinh(bl(j)*x/L),3);
end
%%%%%%%%%%total deflection
z=0;
for j=1:10
    z=z+qm(j)*phim(j);
end
tv=linspace(0,20,20);
z=matlabFunction(z);
for i=1:length(xv)
    for j=1:length(tv)
        zplot(i,j)=z(tv(j),xv(i));
    end
end

figure
mesh(real(zplot))

%%%%%%%%%% shear stress
y=h/2;
sig=0;
for j=1:10
    sig=sig+qm(j)*phid(j);
end
sig=matlabFunction(sig*(E*I)/(I/y));

for i=1:length(xv)
    for j=1:length(tv)
        sigplot(i,j)=sig(tv(j),xv(i));
    end
end
figure
mesh(real(sigplot))
%%%%%%%%%%bending moment
figure
mesh(real(sigplot*(I/y)))
%%%%%%%%%%bending stress
Q=B*h^2/2;%change
tau=0;
for j=1:10
    tau=tau+qm(j)*phi3d(j);
end
tau=matlabFunction(tau*(E*I)/(I*B/Q));

for i=1:length(xv)
    for j=1:length(tv)
        tauplot(i,j)=tau(tv(j),xv(i));
    end
end
figure
mesh(real(tauplot))
toc