  clear all;  clc;
  
%% inputs
global mu 
file='As1.mat';       % As1.mat/mt1.mat/p6.mat/ri1.mat/rs2_1.mat(high EL)/rs2_2.mat
k=100; q=600;           % kth and qth vectors used in Gibbs determination
load(file)              % load the corresponding mat file
T=0;
deg = pi/180;
mu = 398600;
j=size(A);                        % number of readings
%% Time, month, date, true anomaly
ut=hh+mm/60+ss/60/60+mss/1000/60/60;    % calculate UT hours
for i=1:length(Year)
    [m(i),d(i)] = md(Year(i),Day(i));
end
m=m';d=d';                              % month and date
theta=LST(Year, m, d, ut,EL);           % local sidereal time

%% Lambert's method

 
r1=r_from_observe(rho(k), A(k), a(k), theta(k), phi, H);         % initial postition vector
r2=r_from_observe(rho(q), A(q), a(q), theta(q), phi, H);         % final position vector
dt=3600*(ut(q)-ut(k));                       % time interval b/w final and initial

string='pro';
[v1,v2] = lambert(r1, r2, dt, string);
coe = coe_from_sv(r1, v1, mu); 
TA1 = coe(6);
if coe(7)<0
    string='retro';
    [v1,v2] = lambert(r1, r2, dt, string);
    coe = coe_from_sv(r1, v1, mu);
end
coe = coe_from_sv(r2, v2, mu);                
TA2 = coe(6);

%% Print outputs

fprintf('-----------------------------------------------------')
fprintf('\n\n Input data:\n');
fprintf('\n Gravitational parameter (km^3/s^2)  = %g\n', mu);
fprintf('\n r1 (km)                             = [%g %g %g]', ...
r1(1), r1(2), r1(3))
fprintf('\n r2 (km)                             = [%g %g %g]', ...
r2(1), r2(2), r2(3))
fprintf('\n Elapsed time (s)                    = %g', dt);
fprintf('\n\n Lambert OD Solution:\n')
fprintf('\n v1 (km/s)                           = [%g %g %g]', ...
v1(1), v1(2), v1(3))
fprintf('\n v2 (km/s)                           = [%g %g %g]', ...
v2(1), v2(2), v2(3))
fprintf('\n\n Orbital elements:')
fprintf('\n Angular momentum (km^2/s)           = %g', coe(1))
fprintf('\n Eccentricity                        = %g', coe(2))
fprintf('\n Inclination (deg)                   = %g', coe(4)/deg)
fprintf('\n Type of orbit Prograde / Retrograde = %c%c%c%c%c', string)
fprintf('\n RA of ascending node (deg)          = %g', coe(3)/deg)
fprintf('\n Argument of perigee (deg)           = %g', coe(5)/deg)
fprintf('\n True anomaly initial (deg)          = %g', TA1/deg)
fprintf('\n True anomaly final (deg)            = %g', TA2/deg)
fprintf('\n Semimajor axis (km)                 = %g', coe(7))
fprintf('\n Periapse radius (km)                = %g', coe(1)^2/mu/(1 +coe(2)))
%...If the orbit is an ellipse, output its period:
if coe(2)<1
T = 2*pi/sqrt(mu)*coe(7)^1.5;
fprintf('\n Period                              = %g minutes', T/60)
end
fprintf('\n-----------------------------------------------------\n')
%% Orbit  Plot
m1 = 5.974e24 ; 
R = 6378.178;
m2 = 1000;
t0 = 0;
tf = T/60;  % hours
orbit(m1,m2,R,r2,v2,t0,tf)



