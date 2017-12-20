function [ X,Y,Z,T ] = nucleardecay( N0,a_X,a_Y,t,dt)
%given Y_init=0 and z_init=0

%disceretisation of time from initial to final
T = 0:dt:t;

%arrays to store values of isotopes at each timestep
X = zeros(1,length(T));
Y = zeros(1,length(T));
Z = zeros(1,length(T));

%inital value of x 
X(1) = N0;

%probalities of decay of X and Y
P_X = 1 - exp(-a_X*dt);
P_Y = 1 - exp(-a_Y*dt);

%calculation of isotopevalues at each time step
for i=1:length(T)
    if i~=1
        X(i) = X(i-1) * exp(-a_X*dt);
        Y(i) = Y(i-1) + X(i-1) * P_X - Y(i-1) * P_Y;
        Z(i) = Z(i-1) + Y(i-1) * P_Y;
    end
end
%animations
lx=animatedline( 'Color','red');
axis([0 T(end) 0 N0])
ly=animatedline('Color','blue');
axis([0 T(end) 0 N0])
lz=animatedline('Color','green');
axis([0 T(end) 0 N0])

for i=1:length(T)
    addpoints(lx,T(i),X(i))
    
    addpoints(ly,T(i),Y(i))
    
    addpoints(lz,T(i),Z(i))
    
    pause(0.5)
end


end

