function [  ] = lissajousplot( a,b,d)
%takes input a,b,delta and gives the lissajous plot
%defining parametric equations
syms  t
x=@(t) sin(a*t);
y=@(t) sin(b*t+d);
%for finding which function has maximum time period for sin(a*t+b) time
%period is 2*pi/a so 'a' should be small for maximum time period
if a>b
    c=b;
else
    c=a;
end
T=linspace(0,5*(2*pi/c),1000);%time steps from 0 to 5 times to that of the maximum time period of the two functions for plotting

%using loops to generate x and y values for the timesteps
for i=1:length(T)
    x1(i)=x(T(i));
    y1(i)=y(T(i));
end
figure
line(x1,y1)
%creating lables and manipulating axes loction for figure
title(strcat(['a=',' ',num2str(a),' ','b=',num2str(b),' ','delta=',num2str(d)]))
xlabel('X')
ylabel('Y')
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
end

