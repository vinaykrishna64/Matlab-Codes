%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu
function [ f ] = force2( t,rt,p0 )

if t<=rt
        f=(p0/(rt))*t;
else
       f=p0;
    
        
    end

end

