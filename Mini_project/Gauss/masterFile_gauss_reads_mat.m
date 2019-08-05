% ~~~~~~~~~~~~~~~~~~~~~-~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Example_5_10
% ~~~~~~~~~~~~
%
% This program uses Algorithms 5.4 and 4.2 to obtain the orbital
% elements from the observational data provided in Example 5.10.
%
% deg    - conversion factor between degrees and radians
% pi     - 3.1415926...
% mu     - gravitational parameter (km^3/s^2)

% Re     - equatorial radius of the earth (km)
% f      - earth's flattening factor
% wE     - angular velocity of the earth (rad/s)
% omega  - earth's angular velocity vector (rad/s) in the
%          geocentric equatorial frame

% rho    - slant range of object (km)
% rhodot - range rate (km/s)
% A      - azimuth (deg) of object relative to observation site
% Adot   - time rate of change of azimuth (deg/s)
% a      - elevation angle (deg) of object relative to observation site
% adot   - time rate of change of elevation angle (degrees/s)

% theta  - local sidereal time (deg) of tracking site
% phi    - geodetic latitude (deg) of site
% H      - elevation of site (km)

% r      - geocentric equatorial position vector of object (km)
% v      - geocentric equatorial velocity vector of object (km)

% coe    - orbital elements [h e RA incl w TA a]
%          where
%              h    = angular momentum (km^2/s)
%              e    = eccentricity
%              RA   = right ascension of the ascending node (rad)
%              incl = inclination of the orbit (rad)
%              w    = argument of perigee (rad)
%              TA   = true anomaly (rad)
%              a    = semimajor axis (km)
% rp    - perigee radius (km)
% T     - period of elliptical orbit (s)
%
% User M-functions required: rv_from_observe, coe_from_sv
% --------------------------------------------------------------------

clear;
clc;
close all;
global  f Re wE mu

deg    = pi/180;
f      = 1/298.256421867;
Re     = 6378.13655;
wE     = 7.292115e-5;
mu     = 398600.4418;
eEarth = sqrt(1-(1-f)^2);
%...Data declaration for Example 5.10:
% rho    = 2551;
% rhodot = 0;
% A      = 90;
% Adot   = 0.1130;
% a      = 30;
% adot   = 0.05651;
% myData = xlsread('as1_alldata.xls');
load('mydataMT.mat');
Days = myData(1,3);
Y = myData(1,2);
if mod(Y, 4) == 0
    status = true;
    if mod(Y, 100) == 0
        status=false;
        if mod(Y,400)==0
            status = true;
        end
    end
else
    status=false;
end
if status==false && Days>59
    Days = Days+1;
end
sumdays = 0;
arr = [31 60 90 120 151 181 212 243 273 304 334 366];
i = 1;
while(Days>arr(i))
    i = i+1;
end
M = i; %neglecting february and 31 days months
D = Days - arr(i-1);
lon  = 77.5109; %BL2 longitude
phi    = deg2rad(13.0344722);  %BL2 latitude
coe = 9000*ones(7,1);
H      = 0.83991; %station height
Nsite = Re/sqrt(1-(eEarth*sin(phi))^2);
tic
while(coe(2) > 0.001 || isnan(coe(2)))
    pts = sort(randperm(size(myData,1),3));
    UT = myData(pts,4)+myData(pts,5)/60+myData(pts,6)/3600+myData(pts,7)/3600000;
    theta = LST(Y, M, D, UT, lon);
    theta = deg2rad(theta);
    rho = myData(pts,9);
    rhodot = myData(pts,23);
    A = myData(pts,21);
    Adot = myData(pts,24);
    a = myData(pts,22);
    adot = myData(pts,25);
    %...Algorithm 5.4:
    for i=1:3
        Rsite(i,:) = [(Nsite+H)*cos(phi)*cos(theta(i)),(Nsite+H)*cos(phi).*sin(theta(i)),(Nsite*(1-eEarth^2)+H)*sin(phi)];
        [rR(:,i),vV(:,i)] = rv_from_observe(rho(i), rhodot(i), A(i), Adot(i), a(i), adot(i), rad2deg(theta(i)), rad2deg(phi), H);
    end
    
    %Check coplanarity
    %Gibbs
    [v2, ierr] = gibbs(rR(:,1), rR(:,2), rR(:,3));
    
    %...If the vectors r1, r2, r3, are not coplanar, abort:
    if ierr == 1
        fprintf('\n  These vectors are not coplanar.\n\n');
        continue;
    end
    
    %Gauss
    for i=1:3
        [ra(i),dec(i)] = ra_and_dec_from_az_el(A(i),a(i),phi/deg,theta(i)/deg);
        ra(i) = deg2rad(ra(i));
        dec(i) = deg2rad(dec(i));
        rhoT(i,:) = [cos(dec(i))*cos(ra(i)); cos(dec(i))*sin(ra(i)); sin(dec(i))];
    end
    [r, v, r_old, v_old] = gauss(rhoT(1,:), rhoT(2,:), rhoT(3,:), ...
        Rsite(1,:),    Rsite(2,:),  Rsite(3,:), ...
        UT(1)*3600,      UT(2)*3600,    UT(3)*3600);
    
    
    %...Algorithm 4.2 for the initial estimate of the state vector
    %   and for the iteratively improved one:
    coe_old = coe_from_sv(r_old,v_old,mu);
    coe     = coe_from_sv(r,v,mu);
    
end
toc
%...Echo the input data and output the solution to
%   the command window:
fprintf('-----------------------------------------------------')
fprintf('\n Example 5.11: Orbit determination by the Gauss method\n')
fprintf('\n Radius of earth (km)               = %g', Re)
fprintf('\n Flattening factor                  = %g', f)
fprintf('\n Gravitational parameter (km^3/s^2) = %g', mu)
fprintf('\n\n Input data:\n');
fprintf('\n Latitude (deg)                = %g', phi/deg);
fprintf('\n Altitude above sea level (km) = %g', H);
fprintf('\n\n Observations:')
fprintf('\n               Right')
fprintf('                                     Local')
fprintf('\n   Time (s)   Ascension (deg)   Declination (deg)')
fprintf('   Sidereal time (deg)')
for i = 1:3
    fprintf('\n %9.4g %11.4f %19.4f %20.4f', ...
        UT(i), ra(i)/deg, dec(i)/deg, theta(i)/deg)
end

fprintf('\n\n Solution:\n')

fprintf('\n Without iterative improvement...\n')
fprintf('\n');
fprintf('\n r (km)                          = [%g, %g, %g]', ...
    r_old(1), r_old(2), r_old(3))
fprintf('\n v (km/s)                        = [%g, %g, %g]', ...
    v_old(1), v_old(2), v_old(3))
fprintf('\n');

fprintf('\n   Angular momentum (km^2/s)     = %g', coe_old(1))
fprintf('\n   Eccentricity                  = %g', coe_old(2))
fprintf('\n   RA of ascending node (deg)    = %g', coe_old(3)/deg)
fprintf('\n   Inclination (deg)             = %g', coe_old(4)/deg)
fprintf('\n   Argument of perigee (deg)     = %g', coe_old(5)/deg)
fprintf('\n   True anomaly (deg)            = %g', coe_old(6)/deg)
fprintf('\n   Semimajor axis (km)           = %g', coe_old(7))
fprintf('\n   Periapse radius (km)          = %g', coe_old(1)^2 ...
    /mu/(1 + coe_old(2)))
%...If the orbit is an ellipse, output the period:
if coe_old(2)<1
    T = 2*pi/sqrt(mu)*coe_old(7)^1.5;
    fprintf('\n   Period:')
    fprintf('\n     Seconds                     = %g', T)
    fprintf('\n     Minutes                     = %g', T/60)
    fprintf('\n     Hours                       = %g', T/3600)
    fprintf('\n     Days                        = %g', T/24/3600)
end

fprintf('\n\n With iterative improvement...\n')
fprintf('\n');
fprintf('\n r (km)                          = [%g, %g, %g]', r(1), r(2), r(3))
fprintf('\n v (km/s)                        = [%g, %g, %g]', v(1), v(2), v(3))
fprintf('\n');
fprintf('\n   Angular momentum (km^2/s)     = %g', coe(1))
fprintf('\n   Eccentricity                  = %g', coe(2))
fprintf('\n   RA of ascending node (deg)    = %g', coe(3)/deg)
fprintf('\n   Inclination (deg)             = %g', coe(4)/deg)
fprintf('\n   Argument of perigee (deg)     = %g', coe(5)/deg)
fprintf('\n   True anomaly (deg)            = %g', coe(6)/deg)
fprintf('\n   Semimajor axis (km)           = %g', coe(7))
fprintf('\n   Periapse radius (km)          = %g', coe(1)^2 ...
    /mu/(1 + coe(2)))
%...If the orbit is an ellipse, output the period:
if coe(2)<1
    T = 2*pi/sqrt(mu)*coe(7)^1.5;
    fprintf('\n   Period:')
    fprintf('\n     Seconds                     = %g', T)
    fprintf('\n     Minutes                     = %g', T/60)
    fprintf('\n     Hours                       = %g', T/3600)
    fprintf('\n     Days                        = %g', T/24/3600)
end
%%
fid=fopen('Gauss.txt','w');
fprintf(fid, '\n Example 5.11: Orbit determination by the Gauss method\n');
fprintf(fid, '\n Radius of earth (km)               = %g', Re);
fprintf(fid, '\n Flattening factor                  = %g', f);
fprintf(fid, '\n Gravitational parameter (km^3/s^2) = %g', mu);
fprintf(fid, '\n\n Input data:\n');
fprintf(fid, '\n Latitude (deg)                = %g', phi/deg);
fprintf(fid, '\n Altitude above sea level (km) = %g', H);
fprintf(fid, '\n\n Observations:');
fprintf(fid, '\n               Right');
fprintf(fid, '                                     Local');
fprintf(fid, '\n   Time (s)   Ascension (deg)   Declination (deg)');
fprintf(fid, '   Sidereal time (deg)');
for i = 1:3
    fprintf(fid, '\n %9.4g %11.4f %19.4f %20.4f', UT(i), ra(i)/deg, dec(i)/deg, theta(i)/deg);
end

fprintf(fid, '\n\n Solution:\n');

fprintf(fid, '\n Without iterative improvement...\n');
fprintf(fid, '\n');
fprintf(fid, '\n r (km)                          = [%g, %g, %g]', r_old(1), r_old(2), r_old(3));
fprintf(fid, '\n v (km/s)                        = [%g, %g, %g]', v_old(1), v_old(2), v_old(3));
fprintf(fid, '\n');

fprintf(fid, '\n   Angular momentum (km^2/s)     = %g', coe_old(1));
fprintf(fid, '\n   Eccentricity                  = %g', coe_old(2));
fprintf(fid, '\n   RA of ascending node (deg)    = %g', coe_old(3)/deg);
fprintf(fid, '\n   Inclination (deg)             = %g', coe_old(4)/deg);
fprintf(fid, '\n   Argument of perigee (deg)     = %g', coe_old(5)/deg);
fprintf(fid, '\n   True anomaly (deg)            = %g', coe_old(6)/deg);
fprintf(fid, '\n   Semimajor axis (km)           = %g', coe_old(7));
fprintf(fid, '\n   Periapse radius (km)          = %g', coe_old(1)^2 /mu/(1 + coe_old(2)));
%...If the orbit is an ellipse, output the period:
if coe_old(2)<1
    T = 2*pi/sqrt(mu)*coe_old(7)^1.5;
    fprintf(fid, '\n   Period:');
    fprintf(fid, '\n     Seconds                     = %g', T);
    fprintf(fid, '\n     Minutes                     = %g', T/60);
    fprintf(fid, '\n     Hours                       = %g', T/3600);
    fprintf(fid, '\n     Days                        = %g', T/24/3600);
end

fprintf(fid, '\n\n With iterative improvement...\n');
fprintf(fid, '\n');
fprintf(fid, '\n r (km)                          = [%g, %g, %g]', r(1), r(2), r(3));
fprintf(fid, '\n v (km/s)                        = [%g, %g, %g]', v(1), v(2), v(3));
fprintf(fid, '\n');
fprintf(fid, '\n   Angular momentum (km^2/s)     = %g', coe(1));
fprintf(fid, '\n   Eccentricity                  = %g', coe(2));
fprintf(fid, '\n   RA of ascending node (deg)    = %g', coe(3)/deg);
fprintf(fid, '\n   Inclination (deg)             = %g', coe(4)/deg);
fprintf(fid, '\n   Argument of perigee (deg)     = %g', coe(5)/deg);
fprintf(fid, '\n   True anomaly (deg)            = %g', coe(6)/deg);
fprintf(fid, '\n   Semimajor axis (km)           = %g', coe(7));
fprintf(fid, '\n   Periapse radius (km)          = %g', coe(1)^2 /mu/(1 + coe(2)));
%...If the orbit is an ellipse, output the period:
if coe(2)<1
    T = 2*pi/sqrt(mu)*coe(7)^1.5;
    fprintf(fid, '\n   Period:');
    fprintf(fid, '\n     Seconds                     = %g', T);
    fprintf(fid, '\n     Minutes                     = %g', T/60);
    fprintf(fid, '\n     Hours                       = %g', T/3600);
    fprintf(fid, '\n     Days                        = %g', T/24/3600);
end
fclose(fid);
%%
fprintf('\n-----------------------------------------------------\n')

save('AS1_GAUSS2.mat');