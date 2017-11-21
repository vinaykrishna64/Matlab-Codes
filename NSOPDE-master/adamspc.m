function [ ] = adamspc( )
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

%initalvalues using RK 4th order 

k1=f(u(1),t(1));
k2=f(u(1)+(h*k1/2),t(1)+h/2);
k3=f(u(1)+(h*k2/2),t(1)+h/2);
k4=f(u(1)+k3,t(1)+h);
u(2)=u(1)+(h*(k1+2*k2+2*k3+k4))/6;

%predictor corrector set

%iteration control
i=3;
disp('specify how to iterate')
ask=input('if you have a bound for difference give y if not give n  ');
 if ask=='y'
  eps=input('enter the bounded value epsilon\');
 else
 count=input('how many time do you want to iterate the corrector\');
 end

while 1
%predictor
u(i)=u(i-1)+(3*h/2)*f(u(i-1),t(i-1))+(-h/2)*f(u(i-2),t(i-2));
%corrector
if ask=='y'
   while 1
     k=u(i)
     u(i)=u(i-1)+h*((3/2)*f(u(i-1),t(i-1))+(-1)*f(u(i-2),t(i-2))+(5/12)*f(u(i),t(i)));
     if  abs(u(i)-k)>eps
        continue
     else
         break
     end
    end
else
    for iter=1:count
        u(i)=u(i-1)+h*((3/2)*f(u(i-1),t(i-1))+(-1)*f(u(i-2),t(i-2))+(5/12)*f(u(i),t(i)));
        
    end
end
i=i+1;

if i>n
    break
else
    continue
end

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

