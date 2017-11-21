%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu
clear all
%given data and initial conditions
m=16;
k=339;
w=sqrt(k/m);
zeta=[0.02, 0.05, 0.10, 0.20, 0.50];
x_0=[1 ];
v_0=[0 ];

for i=1:length(zeta)
    p=zeta(i);%additional as zeta is too big to type
    figure
    for j=1:length(x_0)
        if p<1 %case 1 where zeta <1
           B2= (v_0(j)+p*w*x_0(j))/(w*sqrt(1-p^2));
            syms t
            xf=@(t) exp(-p*w*t)*(x_0(j)*cos(w*t*sqrt(1-p^2))+B2*sin(w*t*sqrt(1-p^2)));
            
        else %case 2 where zeta >=1
            syms t
            freq=-1*(p*w+w*sqrt(p^2-1));
            xf=@(t) exp(freq*t)*x_0(j);
            
          end
        %x interval
        xc=[0,40*2*pi/w];
        
        titlestr=strcat('zeta','=',num2str(zeta(i)),'...','IC',num2str(j));
        %creating tabs
        
        tab1 = uitab('Title','DISPLACEMENT');
        axes(tab1)
        fplot(xf,xc)
        xlabel('time(s)')
        ylabel('displacement(m)')
        title(titlestr)
        
  
    end   
end
