%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  15NA10016 N.S.ViNAY Krishna Rayudu
function [f]=force(t,alp,rt,p0)
    if t<=alp*rt
        f=(p0/(alp*rt))*t;
    elseif t>alp*rt & t<=rt
        f=p0-((0.9*p0/(rt*(1-alp)))*(t-rt*alp));
    else
        f=p0*0.1;
    end
end

