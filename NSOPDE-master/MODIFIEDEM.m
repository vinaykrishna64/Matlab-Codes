function [ ] = MODIFIEDEM( )
%input

disp('*********give RHS of the ODE as a string********* ')
f1=input('f(u,t)  =');
v=input('Initialvalue  =');
a=input('starting point of the interval  =');
b=input('ending point of the interval  =');
h=input('stepsize  =');
eps=input('epsilon = ');
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

for j=2:n
    u(j)=u(j-1)+h*f(u(j-1),t(j-1));%forward euler first approximation
    while 1
        prev=u(j);
         u(j)=u(j-1)+(h/2)*(f(u(j-1),t(j-1))+f(u(j),t(j)));%iteration
        if (abs(u(j)-prev)>eps)%checking modulus
            continue;
        else
            break;
        end
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

