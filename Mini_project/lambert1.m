function  coe=lambert1(r1,r2,dt,string)
global mu
deg = pi/180;
mu = 398600;
[v1, v2] = lambert(r1, r2, dt, string);
coe = coe_from_sv(r1, v1, mu);
TA1 = coe(6);
coe = coe_from_sv(r2, v2, mu);
TA2 = coe(6);
%a=coe(7);
fprintf('-----------------------------------------------------')
fprintf('\n\n Input data:\n');
fprintf('\n Gravitational parameter (km^3/s^2) = %g\n', mu);
fprintf('\n r1 (km) = [%g %g %g]', ...
r1(1), r1(2), r1(3))
fprintf('\n r2 (km) = [%g %g %g]', ...
r2(1), r2(2), r2(3))
fprintf('\n Elapsed time (s) = %g', dt);
fprintf('\n\n Lambert OD Solution:\n')
fprintf('\n v1 (km/s) = [%g %g %g]', ...
v1(1), v1(2), v1(3))
fprintf('\n v2 (km/s) = [%g %g %g]', ...
v2(1), v2(2), v2(3))
fprintf('\n\n Orbital elements:')
fprintf('\n Angular momentum (km^2/s)  = %g', coe(1))
fprintf('\n Eccentricity  = %g', coe(2))
fprintf('\n Inclination (deg) = %g', coe(4)/deg)
fprintf('\n RA of ascending node (deg) = %g', coe(3)/deg)
fprintf('\n Argument of perigee (deg) = %g', coe(5)/deg)
fprintf('\n True anomaly initial (deg) = %g', TA1/deg)
fprintf('\n True anomaly final (deg) = %g', TA2/deg)
fprintf('\n Semimajor axis (km) = %g', coe(7))
fprintf('\n Periapse radius (km) = %g', coe(1)^2/mu/(1 +coe(2)))
%...If the orbit is an ellipse, output its period:
if coe(2)<1
T = 2*pi/sqrt(mu)*coe(7)^1.5;
fprintf('\n Period= %g minutes', T/60)
end
fprintf('\n-----------------------------------------------------\n')