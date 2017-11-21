%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu
clear all
%data

m=16;
k=33;
w=sqrt(k/m);
zeta=0.05;
C=zeta*2*sqrt(k*m);
alp=0.5;
if zeta==0
T=2*pi/w;
else
T=2*pi/(w*sqrt(1-zeta^2));
end
p0=100;


rtv=linspace(0,3*T,100);

for g=1:length(rtv)
vu=[0 ;0];
rt=rtv(g);
dt=rt*0.01;
A=[1/dt k/2; 1/2 (m/dt)+(C/2)];  
det=-k^2/4+(C/(2*dt))+(m/dt^2);
A=A*(1/det);
t=0;
n=1;  
A1=[1/2*dt k/4; 1/4 (m/(2*dt)+C/4)];
while t <= rt*3
    f1=force(t,alp,rt,p0);
    f2=force(t+dt,alp,rt,p0);
%     f1=force2(t,rt,p0);
%     f2=force2(t+dt,rt,p0);
    
    B=A1*[f1 ; f2]*(1/det);
    vu(1:2,n+1)=A*vu(1:2,n)+B;
    n=n+1;
    t=t+dt;
    if dt==0
        break
    end
end
% figure
% plot(linspace(0,t-dt,length(vu(1,1:end))),vu(2,1:end))
if g==1
    dlf(g)=0;
else
dlf(g)=max(vu(2,100:end))/(p0/k);
end
end
figure
plot(rtv/T,dlf)