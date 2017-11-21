function [  ] = RKIM4( )
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
for j=2:n
    x1=t(j-1)+(((3-3^(1/2))/6)*h);
    x2=t(j-1)+(((3+3^(1/2))/6)*h);
    syms k1v k2v
    eq1=k1v==h*f(u(j-1)+k1v/4+((1/4)-(3^.5/6))*k2v,x1);
    eq2=k2v==h*f(u(j-1)+k2v/4+((1/4)+(3^.5/6))*k1v,x2);
    [k1sc,k2sc]=solve([eq1,eq2],[k1v,k2v]); 
    
    %problem exit due to all complex roots
    len=length(k1sc);
    for i=1:len
        tf1(i)=isreal(k1sc(i));
    end
    len=length(k2sc);
    for i=1:len
        tf(i)=isreal(k2sc(i));
    end
    tf1=sort(tf1);
    tf=sort(tf);
    if tf1(1)==0 || tf(1)==0
        error('all complex roots for some k1 or k2 problem not solvable')
    end
    
   %get solutions for k1 and k2
    len=length(k1sc);
    for i=1:len
        tf1(i)=isreal(k1sc(i));
    end
    count=1;
    for i=1:len
        if tf1(i)~=0
            k1s(count)=k1sc(i);
            count=count+1;
        end
    end
    len=length(k1s);
     k1s=double(k1s);
     
    if len>1
    for i=1:len
        a(i)=abs(k1s(i)-u(j-1));
    end
    a=sort(a);
    for i=1:len
        if a(1)==abs(k1s(i)-u(j-1))
            k1=k1s(i);
            break
        else
            continue
        end
    end
    else
        k1=k1s(1);
    end
    len=length(k2sc);
    for i=1:len
        tf(i)=isreal(k2sc(i));
    end
    count=1;
    for i=1:len
        if tf(i)~=0
            k2s(count)=k2sc(i);
            count=count+1;
        end
    end
    len=length(k2s);
    k2s=double(k2s);
   
    if len>1
    for i=1:len
        a(i)=abs(k1s(i)-k1);
    end
    a=sort(a);
    for i=1:len
        if a(1)==abs(k1s(i)-k1)
            k2=k2s(i);
            break
        else
            continue
        end
    end
    else
        k2=k2s(1);
    end
    u(j)=u(j-1)+(.5*(k1+k2));
    
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

