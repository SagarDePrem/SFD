clc
clear all

load('ri1.mat')
mu =398600;
alt=H*1000;
lon=EL*pi/180;
lat=phi*pi/180; deg=pi/180;
% T= ;% Julian date
for i=1:length(Year)
    [m(i),d(i)] = md(Year(i),Day(i));
end
m=m';d=d';                              % month and date
 
T = julian_date(Year, m , d , hh, mm, ss);
 

j=size(A);

e=1; %Just an initialization for the loop
kk=0; %Loop counter
 coe=[0 0 0 0 0 0 0];
while e>0.027 || isnan(coe(2)) 
    kk    = kk+1
    pts   = sort(randperm(j(1),3));
    p=pts(1); q=pts(2); r=pts(3);
    T1=[T(p);T(q);T(r)];
A1=[A(p) A(q) A(r)];
a1=[a(p) a(q) a(r)];
AZI_ELE=[A1;a1];
 RV = laplace_orbit_fit(lat, lon, alt, T1, AZI_ELE)/1000;
 R=RV(1:3);
 V=RV(4:6);
 coe = coe_from_sv(R,V,mu); 
 e=coe(2);
end
 fprintf('\n-----------------------------------------------------\n')
 fprintf('\n\n Orbital elements:')
fprintf('\n Angular momentum (km^2/s)           = %g', coe(1))
fprintf('\n Eccentricity                        = %g', coe(2))
fprintf('\n Inclination (deg)                   = %g', coe(4)/deg)
fprintf('\n RA of ascending node (deg)          = %g', coe(3)/deg)
fprintf('\n Argument of perigee (deg)           = %g', coe(5)/deg)
fprintf('\n Semimajor axis (km)                 = %g', coe(7))
fprintf('\n Periapse radius (km)                = %g', coe(1)^2/mu/(1 +coe(2)))
%...If the orbit is an ellipse, output its period:
if coe(2)<1
T = 2*pi/sqrt(mu)*coe(7)^1.5;
fprintf('\n Period                              = %g minutes', T/60)
end
fprintf('\n-----------------------------------------------------\n')