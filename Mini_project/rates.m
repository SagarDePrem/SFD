% ~~~~~~~~~~~~~~~~~~~~~~~~
function dydt = rates(t,f)
% ~~~~~~~~~~~~~~~~~~~~~~~~
%{
This function calculates the acceleration vector using Equation 2.22
t - time
f - column vector containing the position vector and the
velocity vector at time t
x, y, z - components of the position vector r
r - the magnitude of the the position vector
vx, vy, vz - components of the velocity vector v
ax, ay, az - components of the acceleration vector a
dydt - column vector containing the velocity and acceleration
components
%}
% ------------------------
x = f(1);
y = f(2);
z = f(3);
vx = f(4);
vy = f(5);
vz = f(6);
r = norm([x y z]);
ax = -mu*x/r^3;
ay = -mu*y/r^3;
az = -mu*z/r^3;
dydt = [vx vy vz ax ay az]';
end %rates