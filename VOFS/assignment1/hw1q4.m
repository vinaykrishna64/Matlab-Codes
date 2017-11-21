%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu
clear all
%given
m=16;
k=339;
w=sqrt(k/m);
zeta=[0 0.02 1 1.2];
lamda=[0.5 1 1.5];
p=10^(-2);%%%%%%% change p
mu=m/p;
e=[0.1 0.5 1];

for o=1:length(e)
for k=1:length(lamda)
    la=lamda(k);
    we=la*w;
    figure
    for l=1:length(zeta)
   
        C=zeta(l)*2*sqrt(k*m);
        phase=atan2(C*we,(k-m*we^2));
        x=@(t) mu*e(o)*we^2/sqrt(((k-m*we^2)^2+(w*C)^2))*sin(we*t-phase);
        
        tab1 = uitab('Title',strcat('lamda',num2str(la),'zeta',num2str(zeta(l))));
        axes(tab1)
        fplot(x,[0 5*2*pi/we])
        
        title(strcat('e',num2str(e(o)),'.....','lamda',num2str(la),'zeta',num2str(zeta(l))))
        
        ft=mu*e(o)*we^2/sqrt(((k-m*we^2)^2+(w*C)^2))*sqrt(k^2+(C*we)^2);
        disp(strcat('force ',' ','transmitted','=',num2str(ft)))
        
    end
end
end