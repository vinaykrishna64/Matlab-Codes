%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu
clear all
%%%%%%%%%%%%%%         PART 2            %%%%%%%%%%%%%%%%%
%GIVEN
n=[2 10];
k1=[3 6];
delta=log(0.1);
m=16;
k=339;
w=sqrt(k/m);
x_0=[1 0 -0.5 ];
v_0=[0 -1  1];
for i=1:length(k1)
syms x
p=solve(2*pi*x-delta/(k1(i))*sqrt(1-x^2)==0,x);
p=double(abs(p));
disp(strcat('damping ratio is','..',num2str(p)))

    wd=w*sqrt(1-p^2);
    figure
   
    for j=1:length(x_0)
        if p<1 %case 1 where zeta <1
           B2= (v_0(j)+p*w*x_0(j))/(w*sqrt(1-p^2));
            syms t
            xf=@(t) exp(-p*w*t)*(x_0(j)*cos(w*t*sqrt(1-p^2))+B2*sin(w*t*sqrt(1-p^2)));
            
         elseif p >1 %case 2 where zeta >1
            syms t
            A=(v_0(j)+(p+sqrt(p^2-1))*w*x_0(j))/(2*w*sqrt(p^2-1));
            B=(-v_0(j)-(p-sqrt(p^2-1))*w*x_0(j))/(2*w*sqrt(p^2-1));
            xf= A*exp((-p+sqrt(p^2-1))*w*t)+B*exp((-p-sqrt(p^2-1))*w*t);
           
            xf=matlabFunction(xf);
           
        else
            xf=@(t) exp(-w*t)*(x_0(j)+(v_0(j)+x_0(j)*w)*t);
            
        end
        
        %x interval
        xc=[0,10*2*pi/w];
        
        titlestr=strcat('zeta','=',num2str(p),'...','IC',num2str(j));
        %creating tabs
        
        tab1 = uitab('Title','DISPLACEMENT');
        axes(tab1)
        fplot(xf,xc)
        xlabel('time(s)')
        ylabel('displacement(m)')
        title(titlestr)
     
        
    end
   
end


