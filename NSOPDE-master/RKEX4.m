function [  ] = RKEX4( )
%input

disp('*********give RHS of the ODE as a string********* ')
f1=input('f(u,t)  =');
v=input('Initialvalue  =');
a=input('starting point of the interval  =');
b=input('ending point of the interval  =');
h=input('stepsize  =');
d=strcat('Du=',f1);
f2=strcat('@(u,t)',f1);
syms u t
f=str2func(f2);
%exactsolution

u0=strcat('u(',num2str(a),')=',num2str(v));
exsol=dsolve(d,u0);
exsol=matlabFunction(exsol);
%nodalpoints

n=((b-a)/h)+1;
t=a:h:b;
u=t;
u(1)=v;

%RK 4th order 
for i=2:n
    k1=f(u(i-1),t(i-1));
k2=f(u(i-1)+(h*k1/2),t(i-1)+h/2);
k3=f(u(i-1)+(h*k2/2),t(i-1)+h/2);
k4=f(u(i-1)+k3,t(i-1)+h);

u(i)=u(i-1)+(h*(k1+2*k2+2*k3+k4))/6;
end
%table creation
u=double(u);
tab=[t;u];
tab=tab';
figure
table=uitable('data',tab,'columnname',{'x','y'})
%plots
figure
x=a:.01:b;
for i=1:((b-a)/.01)+1
    y(i)=exsol(x(i));
end
plot(x,y);
hold on
plot(t,u,':b');
hold on


end

