function [ Ynew] = notaknot( X,Y,Xnew)
%creates a notaknot spline and interpolates for new x values
n1=length(X);
n2=length(Xnew);
%sorting data
[X,index]=sort(X);
Y=Y(index);
Xnew=sort(Xnew);
%finding ci
for i=1:n1
    for j=1:n1
        if i==1
            if j==1
                Ca(i,j)=(X(3)-X(2));
            elseif j==2
                Ca(i,j)=-(X(3)-X(1));
            elseif j==3
                Ca(i,j)=(X(2)-X(1));
            else
                Ca(i,j)=0;
            end
        elseif i==n1
            if j==n1-2
                Ca(i,j)=(X(n1)-X(n1-1));
            elseif j==n1-1
                Ca(i,j)=-(X(n1)-X(n1-2));
            elseif j==n1
                Ca(i,j)=(X(n1-1)-X(n1-2));
            else
                Ca(i,j)=0;
            end
        else
            if j==i-1
                Ca(i,j)=(X(j+1)-X(j));
            elseif j==i
                Ca(i,j)=2*(X(j+1)-X(j-1));
            elseif j==i+1
                Ca(i,j)=(X(j)-X(j-1));
            else
                Ca(i,j)=0;
            end
        end
    end
end
Cb=zeros(n1,1);
for i=1:n1
    if i~=1 & i~=n1
        Cb(i)=3*((Y(i+1)-Y(i))/(X(i+1)-X(i))-(Y(i)-Y(i-1))/(X(i)-X(i-1)));
    end
end
C=linsolve(Ca,Cb);
for i=1:n1-1
    D(i)=(C(i+1)-C(i))/(3*(X(i+1)-X(i)));
end
for i=1:n1-1
    B(i)=((Y(i+1)-Y(i))/(X(i+1)-X(i)))-C(i)*(X(i+1)-X(i))-D(i)*(X(i+1)-X(i))^2;
end
A=Y(1:n1-1);
Ys(1)=Y(1);
figure
for i=2:n1
    p1=fplot(@(x)A(i-1)+B(i-1)*(x-X(i-1))+C(i-1)*(x-X(i-1))^2+D(i-1)*(x-X(i-1))^3,[X(i-1) X(i)],'k');
    hold on
end
for j=1:n2
    for i=1:n1-1
        if Xnew(j)>X(i) & Xnew(j)<X(i+1)
            Ynew(j)=A(i)+B(i)*(Xnew(j)-X(i))+C(i)*(Xnew(j)-X(i))^2+D(i)*(Xnew(j)-X(i))^3;
        end
    end
end
p2=plot(X,Y,'g*');
hold on
p3=plot(Xnew,Ynew,'r*');
legend([p1,p2,p3],'Spline','Y','Y_{new}');
title('not a knot spline')
xlabel('X axis')
ylabel('Y axis')
end

