%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu
clear all
%given
YA=[0.1 1 10];
lamda=[0.5 1 1.5];
m=16;
k=339;
w=sqrt(k/m);
zeta=[0 0.02 1 1.2];

for l1=1:length(lamda)
    la=lamda(l1);
    we=la*w;
    force=@(t) YA(2)*cos(we*t); %change YA
    figure
    for l=1:length(zeta)
        C=zeta(l)*2*sqrt(k*m);
        phase=atan2(m*C*we^3,(C*we)^2+k*(k-m*we^2))
        x=@(t) YA(2)*sqrt((k^2+(C*we)^2)/((k-m*we^2)^2+(w*C)^2))*cos(we*t-phase);%change YA
        
        tab1 = uitab('Title',strcat('lamda',num2str(la),'zeta',num2str(zeta(l))));
        axes(tab1)
        fplot(force,[0 5*2*pi/we],':')
        hold on
        title(strcat('lamda',num2str(la),'zeta',num2str(zeta(l))))
        fplot(x,[0 5*2*pi/we])
        hold off
        legend('excitation','response')
    end
end
        
        
        