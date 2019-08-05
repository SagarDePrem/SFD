% ~~~~~~~~~~~~~
function output
% ~~~~~~~~~~~~~
%{
This function computes the maximum and minimum radii, the times they
occur and and the speed at those times. It prints those results to
the command window and plots the orbit.
r - magnitude of the position vector at the times in t
imax - the component of r with the largest value
rmax - the largest value of r
imin - the component of r with the smallest value
rmin - the smallest value of r
v_at_rmax - the speed where r = rmax
v_at_rmin - the speed where r = rmin
User subfunction required: light_gray
%}
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
fprintf(' %s\n', datestr(now))
fprintf('\n The initial position is [%g, %g, %g] (km).',...
r0(1), r0(2), r0(3))
fprintf('\n Magnitude = %g km\n', norm(r0))
fprintf('\n The initial velocity is [%g, %g, %g] (km/s).',...
v0(1), v0(2), v0(3))
fprintf('\n Magnitude = %g km/s\n', norm(v0))
fprintf('\n Initial time = %g h.\n Final time = %g h.\n',0,tf/hours)
fprintf('\n The minimum altitude is %g km at time = %g h.',...
rmin-R, t(imin)/hours)
fprintf('\n The speed at that point is %g km/s.\n', v_at_rmin)
fprintf('\n The maximum altitude is %g km at time = %g h.',...
rmax-R, t(imax)/hours)
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