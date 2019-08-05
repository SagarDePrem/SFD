% sunsync2.m      April 24, 2008

% mean orbital inclination of sun-synchronous orbits

% j2 + j4 iterative solution

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% global variables

global j2 j4 j22 mm drdt0 req2 req4 slr2 slr4 ecc2

% astrodynamic and utility constants

j2 = 0.00108263;            % first zonal constant (nd)
j4 = 0.0000016109876;       % fourth zonal constant (nd)
j22 = 0.0000015747419;      % j22 gravity constant (nd)
mu = 398600.5;              % Earth gravitational constant (km^3/sec^2)
req = 6378.14;              % Earth equatorial radius (kilometers)

pidiv2 = 0.5 * pi;          % pi/2
rtd = 180 / pi;             % radians to degrees
dtr = pi / 180;             % degrees to radians
rtol = 0.00000001;          % Brent root-finding tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% begin main program

clc; home;

fprintf('              program sunsync2 \n\n');

fprintf(' < sun-synchronous orbits - j2 + j4 solution > \n\n');

fprintf('<1> input mean semimajor axis and eccentricity \n\n');

fprintf('<2> input mean perigee and apogee altitudes \n\n');

selection = input('selection (1 or 2) ? ');

if (selection == 1)

   fprintf('\n');

   sma = input('please input the mean semimajor axis (kilometers)? ');

   fprintf('\n');

   ecc = input('please input the mean orbital eccentricity (non-dimensional)? ');

   % compute perigee and apogee altitudes

   hp = sma * (1 - ecc) - req;

   ha = sma * (1 + ecc) - req;

elseif (selection == 2)

   fprintf('\n');

   hp = input('please input the mean perigee altitude (kilometers)? ');

   fprintf('\n');

   ha = input('please input the mean apogee altitude (kilometers)?  ');

   % compute perigee and apogee geocentric radii

   rp = req + hp;
   
   ra = req + ha;

   % compute semimajor axis and eccentricity

   sma = 0.5 * (rp + ra);
   
   ecc = (ra - rp) / (ra + rp);
   
end

% required nodal regression rate (radians/second)

drdt0 = dtr * (360 / 365.2422) / 86400;

% Keplerian mean motion (radians/second)

mm = sqrt(mu / (sma * sma * sma));

% orbital semiparameter (kilometers)

slr = sma * (1 - ecc * ecc);
   
% inclination initial guess

ipar0 = -(2 / 3) * slr * slr * drdt0 / (j2 * req * req * mm);

% check for possible solution

if (abs(ipar0) > 1)
    
   % no solution
     
   clc; home;
   
   disp ('no sun-synchronous solution !!');
   
   pause;
   
   xroot = 0;
   
else
    
   % calculate solution

   slr2 = slr * slr;
   slr4 = slr2 * slr2;

   req2 = req * req;
   req4 = req2 * req2;

   ecc2 = ecc * ecc;

   % inclination initial guess

   inc0 = acos(ipar0);

   % define inclination bracket

   inc1 = inc0 - 1 * dtr;
    
   if (inc1 < pidiv2) 
      inc1 = pidiv2;
   end

   inc2 = inc0 + 1 * dtr;

   if (inc2 > pi)
      inc2 = pi;
   end

   % solve nonlinear inclination equation

   [xroot, froot] = brent('ss2func', inc1, inc2, rtol);

end

% print results

clc; home;

fprintf('                program sunsync2 \n\n');
   
fprintf('   < sun-synchronous orbits - j2 + j4 solution > \n\n');
     
fprintf('mean perigee altitude      %12.4f  kilometers \n\n', hp);
 
fprintf('mean apogee altitude       %12.4f  kilometers \n\n', ha);

fprintf('mean semimajor axis        %12.4f  kilometers \n\n', sma);
     
fprintf('mean orbital eccentricity  %12.10f \n\n', ecc);

fprintf('mean orbital inclination   %12.4f  degrees \n\n', xroot * rtd);




