%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu
clear all
%given data and initial conditions
m=16;
k=339;
w=sqrt(k/m);
zeta=[0 0.02 1 1.2];
x_0=[1 0 -0.5 ];
v_0=[0 -1  1];
lamda=[0.5 1 1.5];
F0=1;
for k=1:length(lamda)
        la=lamda(k);
    for i=1:length(zeta)
        %complmentary function
        p=zeta(i);%additional as zeta is too big to type
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
            %particular integral
            ci=sqrt(-1);
            we=la*w;
            c=2*p*sqrt(k*m);
            xpa= F0/sqrt((k-m*we^2)^2+(we*c)^2);
            alpha=atan2(-we*c,k-m*we^2);
            xpf=@(t) exp(ci*(we)*t+ci*alpha)*xpa;
            vpfexp=diff(exp(ci*(we)*t+ci*alpha)*xpa);
            vpf=matlabFunction(vpfexp);
            apfexp=diff(exp(ci*(we)*t+ci*alpha)*xpa,2);
            apf=matlabFunction(apfexp);
            %force
            force=@(t) F0*exp(ci*we*t);
            %total functions
            xt=@(t) xf(t)+xpf(t);
            vt=@(t) vf(t)+vpf(t);
            at=@(t) af(t)+apf(t);
            %x interval
            xc=[0,10*2*pi/w];
            xc1=linspace(0,10*2*pi/w,2000);
            for l=1:length(xc1)
                xtv(l)=xt(xc1(l));
                vtv(l)=vt(xc1(l));
                atv(l)=at(xc1(l));
                forcev(l)=force(xc1(l));
            end
            titlestr=strcat('lambda',num2str(la),'..','zeta','=',num2str(zeta(i)),'...','IC',num2str(j));
            %creating tabs

            tab1 = uitab('Title','DISPLACEMENT');
            axes(tab1)
            plot(xc1,xtv,'b',xc1,forcev,'-.k')
            xlabel('time(s)')
            ylabel('displacement(m)')
            title(titlestr)
            tab2 = uitab('Title','VELOCITY');
            axes(tab2)
            plot(xc1,vtv,'b',xc1,forcev,'-.k')
            xlabel('time(s)')
            ylabel('velocity(m/s)')
            title(titlestr)
            tab3 = uitab('Title','ACCELERATION');
            axes(tab3)
            plot(xc1,atv,'b',xc1,forcev,'-.k')
            xlabel('time(s)')
            ylabel('acceleration(m/s^2)')
            title(titlestr)
            
        end
   
    end
end
