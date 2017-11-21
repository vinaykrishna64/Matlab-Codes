
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
