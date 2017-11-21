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


rtv=linspace(0,20*T,10);

for g=1:length(rtv)
vu=[0 ;0];
rt=rtv(g);
dt=rt*0.01;
A=[1 -k*dt/m; dt 1+C*dt/m];  
det=1+(C*dt/m)+(k*dt^2/m);
A=A*(1/det);
t=0;
n=1;  

while t <= rt*3
%     f1=force(t+dt,alp,rt,p0);
    f1=force2(t,rt,p0);
    
    B=A*[dt*f1/m ;0];
    vu(1:2,n+1)=A*vu(1:2,n)+B;
    n=n+1;
    t=t+dt;
    if dt==0
        break
    end
end
figure
plot(linspace(0,t-dt,length(vu(1,1:end))),vu(2,1:end))
if g==1
    dlf(g)=0;
else
dlf(g)=max(vu(2,100:end))/(p0/k);
end
end
figure
plot(rtv/T,dlf)