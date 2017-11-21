
function []=BEM()
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

%y_i+1 = y_i + h * f_i+1
for j=2:n
    syms y
    eq= y==double(u(j-1))+h*f(y,t(j));
    
    yi=double(solve(eq,y));%gives many values
    count=1;
    for i=1:length(yi)
        if isreal(yi(i)) ~=0
            y(count)=yi(i);
            count=count+1;
        end
    end
    len=length(y);
    if len>1
    for i=1:len
        a(i)=abs(y(i)-u(j-1));
    end
    a=sort(a);
    for i=1:len
        if a(1)==abs(y(i)-u(j-1))
            u(j)=y(i);
            break
        else
            continue
        end
    end
    else
        u(j)=y(1);
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