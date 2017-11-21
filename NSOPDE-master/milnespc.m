function [ ] = milnespc( )
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
y=t;
y(1)=v;
x=t;
%initalvalues using taylor series order 3 4STEP METHOD
syms u t
g=sym(f);
g2=diff(g,t);
g2=strcat('@(u,t)',char(g2));
g2=str2func(g2);

for w=1:3
    uf=y(w);
    tf=x(w);
y(w+1)=uf+h*f(uf,tf)+(h^2/2)*g2(uf,tf);
end


%predictor corrector set

%iteration control
i=5;
disp('specify how to iterate')
ask=input('if you have a bound for difference give y if not give n  ');
 if ask=='y'
  eps=input('enter the bounded value epsilon\');
 else
 count=input('how many time do you want to iterate the corrector\');
 end

while 1
%predictor
y(i)=y(i-4)+(8*h/3)*f(y(i-3),x(i-3))+(-4*h/3)*f(y(i-2),x(i-2))+(8*h/3)*f(y(i-1),x(i-1));
%corrector
if ask=='y'
   while 1
     k=y(i)
     y(i)=y(i-2)+(h/3)*(4*f(y(i-1),x(i-1))+f(y(i-2),x(i-2))+f(y(i),x(i)));
     if  abs(y(i)-k)>eps
        continue
     else
         break
     end
    end
else
    for iter=1:count
        y(i)=y(i-2)+(h/3)*(4*f(y(i-1),x(i-1))+f(y(i-2),x(i-2))+f(y(i),x(i)));
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
y=double(y);
tab=[x;y];
tab=tab';
figure
table=uitable('data',tab,'columnname',{'x','y'})
%plots
figure
plot(x,y,':b');
hold on
p=a:.01:b;
for i=1:((b-a)/.01)+1
    o(i)=exsol(p(i));
end
plot(p,o);
hold on

end

