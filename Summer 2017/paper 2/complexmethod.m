function [XF ] = complexmethod( A,B )
if iscolumn(B)==0
    B=B';
end
AR=real(A);
AC=-i*(A-real(A));
BR=real(B);
BC=-i*(B-real(B));
%singular value decomposition
% %     [U,S,V]=svd(A,0);
% %     X=V*((U'*B)./diag(S));
n=length(A);
AF=[[AR -AC];[AC AR]];
BF=[BR;BC];
X=AF\BF;
XF=X(1:n)+i*X(n+1:end);
end

