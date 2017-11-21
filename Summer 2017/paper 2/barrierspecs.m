function [ i,T,X,h ] = barrierspecs(  )
disp('Type G / B for gap or bar respectively followed by the interval and Type END when you are done entering intervals  ***all G B END are strings ')
i=1;
while true
    Temp=input('G / B') ;
    if Temp=='END'
        break
    else
    T(i)=Temp;   
    X{i}=input('interval');
    i=i+1;
    end
end
i=i-1;
h=X{i}(2);
end

