function dfdt = rates(t,f)
 
% This function calculates the spacecraft acceleration from its
% position and velocity at time t.
R = f(1:3)';                        %Position vector (km/s)
r =norm(R);                         %Distance from earth’s center (km)
alt = r - RE;                       %Altitude (km)
rho = atmosphere(alt);              %Air density from US Standard Model (kf/m^3)
V = f(4:6)';                        %Velocity vector (km/s)
Vrel = V - cross(wE,R);             %Velocity relative to the atmosphere (km/s)
vrel = norm(Vrel);                  %Speed relative to the atmosphere (km/s)
uv = Vrel/vrel;                     %Relative velocity unit vector
ap = -CD*A/m*rho*...                %Acceleration due to drag (m/s^2)
(1000*vrel)^2/2*uv;                 %(converting units of vrel from km/s to m/s)
a0 = -mu*R/r^3;                     %Gravitational acceleration (km/s^2)
a = a0 + ap/1000;                   %Total acceleration (km/s^2)
dfdt = [V a]';                      %Velocity and the acceleraion returned to ode45
end                                 %rates