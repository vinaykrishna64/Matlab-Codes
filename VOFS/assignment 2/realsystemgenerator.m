% genrate a randomn set of inputs so that you can have a real system to
% work on
%you're welcome :P 15NA10016
ev=i*ones(1,3);  %just not zero
syms wsq
while imag(ev(1))~=0  | imag(ev(2))~=0  | imag(ev(3))~=0    
m=rand([1 3]);
k=rand([1 6]);
M=double([m(1) 0 0;0 m(2) 0 ;0 0 m(3)]);
K=double([k(1)+k(2)+k(5) -k(2) -k(5);-k(2) k(2)+k(3)+k(6) -k(3); -k(5) -k(3) k(3)+k(4)+k(5)]);
%eigen values
ev=solve(det(K-M*wsq)==0);
ev=sort(double(ev));
end
disp('masses')
disp(m)
disp('stifnesses')
disp(k)
disp('eigen values')
disp(ev)
mat=K-wsq*M;
syms b c
matn=mat*[1;b;c];
s=solve([matn(1)==0,matn(2)==0],[b c]);
a2a1=matlabFunction(s.b);
a3a1=matlabFunction(s.c);
for l=1:3
    p(1:3,l)=[1;a2a1(ev(l));a3a1(ev(l))];
end
disp(p)