
clear all
clc
syms wsq
m=[0.6035    0.5261    0.7297];
k=[0.7073    0.7814    0.2880    0.6925    0.5567    0.3965];

M=double([m(1) 0 0;0 m(2) 0 ;0 0 m(3)]);
K=double([k(1)+k(2)+k(5) -k(2) -k(5);-k(2) k(2)+k(3)+k(6) -k(3); -k(5) -k(3) k(3)+k(4)+k(5)]);
%eigen values
ev=solve(det(K-M*wsq)==0);
ev=sort(double(ev));

%frequencies
for l=1:3
    fr(l)=sqrt(ev(l));
    disp(strcat('frequency',num2str(l),' ','=',' ',num2str(double(fr(l)))));
end
%eigen vectors

mat=K-wsq*M;
syms b c
matn=mat*[1;b;c];
s=solve([matn(1)==0,matn(2)==0],[b c]);
a2a1=matlabFunction(s.b);
a3a1=matlabFunction(s.c);
for l=1:3
    p(1:3,l)=[1;a2a1(ev(l));a3a1(ev(l))];
end

% %ptilda
% for l=1:3
%     mod=double(sqrt(p(1:3,l)'*p(1:3,l)));
%     pt(1:3,l)=p(1:3,l)./mod;
% end
GM=double(p'*M*p);
GK=double(p'*K*p);


%%%%%%%%%%%%%               force part                %%%%%%%%%%%%%%%%%%%%
syms f w t
F=[f*cos(w*t);f*sin(w*t);0];
GF=p'*F;

F1=matlabFunction(GF(1));
F2=matlabFunction(GF(2));
fv=1;
wv=10;
interval=5*2*pi/wv;
fplot(@(x)F1(fv,x,wv),[0 interval])
figure
fplot(@(x)F2(fv,x,wv),[0 interval])


% %%%%%%%%%%%%%%%%%     principal co-ordinates  %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ode=strcat(num2str(GK(1,1)),'*D2x','+',num2str(GM(1,1)),'*x','==',char(GF(1)));
% q1=dsolve(ode);
% q1=matlabFunction(q1);
% ode=strcat(num2str(GK(2,2)),'*D2x','+',num2str(GM(2,2)),'*x','==',char(GF(2)));
% q2=dsolve(ode);
% q2=matlabFunction(q2);
% ode=strcat(num2str(GK(3,3)),'*D2x','+',num2str(GM(3,3)),'*x','==',char(GF(3)));
% q3=dsolve(ode);
% q3=matlabFunction(q3);
% syms c1 c2 c3 c4 c5 c6
% q1=q1(c1,c2,fv,t,wv);
% q2=q2(c3,c4,fv,t,wv);
% q3=q3(c5,c6,fv,t,wv);
% 
% %%%%%%%%%%%%%%%%%%%%%%   modal displacements %%%%%%%%%
% 
% xd=p*{[q1];[q2];[q3]};
% 
% xd1=matlabFunction(xd(1));
% xd2=matlabFunction(xd(2));
% xd3=matlabFunction(xd(3));


%%%%%%%%%%% principal co-ordinates %%%%%%%%%    
q=inv(GK-w^2*GM)*GF;
for i=1:3
 q1(i)={matlabFunction(q(i))};
end

for i=1:3
    figure
fplot(@(t) q1{i}(2,t,3),[0 16]);
end

%%%%%%%%%%%%%%%%%%  modal displacements %%%%%%%%%%%
x=p*q;
for i=1:3
 q2(i)={matlabFunction(x(i))};
end


for i=1:3
    figure
fplot(@(t) q2{i}(2,t,3),[0 16]);
end
