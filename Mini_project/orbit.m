function orbit=orbit(m1,m2,R,r0,v0,t0,tf)
%clc; close all; clear all
hours = 3600;
G = 6.6742e-20;
% t0 = 0;
tf=tf*3600;
%...End input data

%...Numerical integration:
mu = G*(m1 + m2);
y0 = [r0 v0]';
[t,y] = rkf45(@rates, [t0 tf], y0) ;

%...Output the results:
output
return
% ~~~~~~~~~~~~~~~~~~~~~~~~
function dydt = rates(t,f)
% ~~~~~~~~~~~~~~~~~~~~~~~~

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
% ~~~~~~~~~~~~~
function output
% ~~~~~~~~~~~~~

% -------------
for i = 1:length(t)
r(i) = norm([y(i,1) y(i,2) y(i,3)]);
end
[rmax imax] = max(r);
[rmin imin] = min(r);
v_at_rmax = norm([y(imax,4) y(imax,5) y(imax,6)]);
v_at_rmin = norm([y(imin,4) y(imin,5) y(imin,6)]);
%...Output to the command window:
fprintf('\n\n--------------------------------------------------------\n')
fprintf('\n Earth Orbit\n')
fprintf('\n The initial position is [%g, %g, %g] (km).',...
r0(1), r0(2), r0(3))
fprintf('\n Magnitude = %g km\n', norm(r0))
fprintf('\n The initial velocity is [%g, %g, %g] (km/s).',...
v0(1), v0(2), v0(3))
fprintf('\n Magnitude = %g km/s\n', norm(v0))
fprintf('\n Initial time = %g min.\n Final time = %g min.\n',0,tf/hours)
fprintf('\n The minimum altitude is %g km at time = %g h.',...
rmin-R, t(imin)/hours/60)
fprintf('\n The speed at that point is %g km/s.\n', v_at_rmin)
fprintf('\n The maximum altitude is %g km at time = %g h.',...
rmax-R, t(imax)/hours/60)
fprintf('\n The speed at that point is %g km/s\n', v_at_rmax)
fprintf('\n--------------------------------------------------------\n\n')
%...Plot the results:
% Draw the planet
[xx, yy, zz] = sphere(100);
surf(R*xx, R*yy, R*zz)
colormap(light_gray)
caxis([-R/100 R/100])
shading interp
% Draw and label the X, Y and Z axes
line([0 2*R], [0 0], [0 0]); text(2*R, 0, 0, 'X')
line( [0 0], [0 2*R], [0 0]); text( 0, 2*R, 0, 'Y')
line( [0 0], [0 0], [0 2*R]); text( 0, 0, 2*R, 'Z')
% Plot the orbit, draw a radial to the starting point
% and label the starting point (o) and the final point (f)
hold on
plot3( y(:,1), y(:,2), y(:,3),'k')
line([0 r0(1)], [0 r0(2)], [0 r0(3)])
text( y(1,1), y(1,2), y(1,3), 'o')
text( y(end,1), y(end,2), y(end,3), 'f')
% Select a view direction (a vector directed outward from the origin)
view([1,1,.4])
% Specify some properties of the graph
grid on
axis equal
xlabel('km')
ylabel('km')
zlabel('km')
% ~~~~~~~~~~~~~~~~~~~~~~~
function map = light_gray
% ~~~~~~~~~~~~~~~~~~~~~~~

% -----------------------
r = 0.8; g = r; b = r;
map = [r g b
0 0 0
r g b];
end %light_gray
end %output
end %orbit
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~