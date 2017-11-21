function [ temp] = gaussseidal(A,B)
N=length(A);
if iscolumn(B)==0
    B=B';
end
if rcond(A)<=0.1
disp('matrix is ill-conditioned result maybe inacurate')
end
temp=zeros(N,1);
%gauss seidel with previous iteration as temp and current as X
temp2=999999;
%x=tx+b form
for p=1:N
    for q=1:N
        if p==q
            T(p,q)=0;
        else 
            T(p,q)=-A(p,q)/A(p,p);
        end
    end
end
iteration=0;
while abs(temp2(1)-temp(1))>=0.0001
    temp2=temp;
    if iscolumn(temp)==0
        temp=temp';
    end
    for p=1:N
        temp(p)=T(p,1:end)*temp+(B(p)/A(p,p));
  
    end
    iteration=iteration+1;
    if iteration>1000
        disp('iteration tolerance of  1000 iterations exceeded ')
        break
    end
end
end

