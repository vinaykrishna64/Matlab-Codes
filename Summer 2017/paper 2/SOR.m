function [ temp ] = SOR( A,B)
if iscolumn(B)==0
    B=B';
end
n=length(A);
L=zeros(n);
D=zeros(n);
U=zeros(n);
%making L U D
for p=1:n
    for q=1:n
        if p<q
            U(p,q)=A(p,q);
        elseif p==q
            D(p,q)=A(p,q);
        else
            L(p,q)=A(p,q);
        end
    end
end
%%%%%%%%    OPTIMAL VALUE OF w B=D-1*(L+U)(JACOBI)
%%%%%%%%    W0PT=2/1+SQRT(1-RHO(b)^2) RHO(b) IS SPECTRAL RADIUS
Bj=inv(D)*(L+U);
RHO=max(abs(eig(Bj)));
if RHO<1
w=2/(1+sqrt(1-RHO^2));
else
    error('SPECTRAL RADIUS OF ITERATION MATRIX>1') 
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp=zeros(n,1);
if iscolumn(temp)==0
   temp=temp';
end
%gauss seidel with previous iteration as temp and current as X
temp2=999999;
%x=(D+wL)^-1 (wB-[wU+(w-1)D]x) form
C1=inv(D+w*L);
C2=(w*U+(w-1)*D);
while abs(temp2(1)-temp(1))>=0.0001
    temp2=temp;
    temp=C1*(w*B-C2*temp);
end
end

