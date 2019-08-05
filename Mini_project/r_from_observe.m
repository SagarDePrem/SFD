function r = r_from_observe(rho, A, a,theta, phi, H)

global f Re wE
f=0.003353; Re=6378.137; wE=72.92*10^-6;
deg = pi/180;
omega = [0 0 wE];
A = A *deg;
a = a *deg;
theta = theta*deg;
phi = phi *deg;
R = [(Re/sqrt(1-(2*f - f*f)*sin(phi)^2) + H)*cos(phi)*cos(theta), (Re/sqrt(1-(2*f - f*f)*sin(phi)^2) + H)*cos(phi)*sin(theta), (Re*(1 - f)^2/sqrt(1-(2*f - f*f)*sin(phi)^2) + H)*sin(phi)];
Rdot = cross(omega, R);
dec = asin(cos(phi)*cos(A)*cos(a) + sin(phi)*sin(a));  
h = acos((cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a))/cos(dec));
if (A > 0) & (A < pi)
h = 2*pi - h;
end
RA = theta - h;
Rho = [cos(RA)*cos(dec) sin(RA)*cos(dec) sin(dec)];
r = R + rho*Rho;
