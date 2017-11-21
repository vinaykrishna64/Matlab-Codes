%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu

clear all
%given data and initial conditions
m=16;
k=339;
w=sqrt(k/m);
zeta=[0 0.02 1 1.2];
x_0=[1 0 -0.5 ];
v_0=[0 -1  1];

for i=1:length(zeta)
    p=zeta(i);%additional as zeta is too big to type
    wd=w*sqrt(1-p^2);
    figure
    for j=1:length(x_0)
        if p<1 %case 1 where zeta <1
           B2= (v_0(j)+p*w*x_0(j))/(w*sqrt(1-p^2));
            syms t
            xf=@(t) exp(-p*w*t)*(x_0(j)*cos(w*t*sqrt(1-p^2))+B2*sin(w*t*sqrt(1-p^2)));
            vfexp= diff(exp(-p*w*t)*(x_0(j)*cos(w*t*sqrt(1-p^2))+B2*sin(w*t*sqrt(1-p^2))));
            vf=matlabFunction(vfexp);
            afexp= diff(exp(-p*w*t)*(x_0(j)*cos(w*t*sqrt(1-p^2))+B2*sin(w*t*sqrt(1-p^2))),2);
            af=matlabFunction(afexp);
        elseif p>1 %case 2 where zeta >1
            syms t
            A=(v_0(j)+(p+sqrt(p^2-1))*w*x_0(j))/(2*w*sqrt(p^2-1));
            B=(-v_0(j)-(p-sqrt(p^2-1))*w*x_0(j))/(2*w*sqrt(p^2-1));
            xf= A*exp((-p+sqrt(p^2-1))*w*t)+B*exp((-p-sqrt(p^2-1))*w*t);
            vfexp= diff(xf);
            vf=matlabFunction(vfexp);
            afexp= diff(vfexp);
            af=matlabFunction(afexp);
            xf=matlabFunction(xf);
            
        else
            xf=@(t) exp(-w*t)*(x_0(j)+(v_0(j)+x_0(j)*w)*t);
            syms t
            vfexp= diff(exp(-w*t)*(x_0(j)+(v_0(j)+x_0(j)*w)*t));
            vf=matlabFunction(vfexp);
            afexp= diff(vfexp);
            af=matlabFunction(afexp);
        end
        
        %x interval
        xc=[0,10*2*pi/w];
        
        titlestr=strcat('zeta','=',num2str(zeta(i)),'...','IC',num2str(j));
        %creating tabs
        
        tab1 = uitab('Title','DISPLACEMENT');
        axes(tab1)
        fplot(xf,xc)
        xlabel('time(s)')
        ylabel('displacement(m)')
        title(titlestr)
        tab2 = uitab('Title','VELOCITY');
        axes(tab2)
        fplot(vf,xc)
        xlabel('time(s)')
        ylabel('velocity(m/s)')
        title(titlestr)
        tab3 = uitab('Title','ACCELERATION');
        axes(tab3)
        fplot(af,xc)
        xlabel('time(s)')
        ylabel('acceleration(m/s^2)')
        title(titlestr)
        if  p<1
            wphase=wd;
        else 
            wphase=w;
        end
        %phasee plane
        %generating data for plots
        xc1=linspace(0,2*2*pi/w,200);
        for l=1:length(xc1)
            xfv(l)=xf(xc1(l));
            vfv(l)=vf(xc1(l));
        end
        tab4 = uitab('Title','PHASE PLANE');    
        axes(tab4)
        plot(xfv,(vfv*(1/wphase)))
        xlabel('displacement(m)')
        ylabel('velocity/w (m)')
        title(titlestr)
    end
   
end
