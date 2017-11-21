function []=F(h,C,a,D,G)
pi=3.1416;
g=9.8;
i=(-1)^(0.5);
N=3;%change to 9
for k=1:length(C)
    W=(g/(C(k)*h))^(0.5);
    kt=roots(C(k),h);
    A=matriceA(kt,W,G,D,a,h);
    B=matriceB(kt,W,G,D,a,h);
   
    X=linsolve(A,B);
    R0(k)=X(1);
end
R0
end

