1% sunsync1.m      April 24, 2008

% mean orbital inclination
% of sun-synchronous orbits

% Kozai j2 perturbations

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% astrodynamic and utility constants

j2 = 0.00108263;             % zonal coefficient (nd)
mu = 398600.5;               % gravitational constant (km^3/sec^2)
req = 6378.14;               % Earth equatorial radius (kilometers)
rtd = 180 / pi;              % radians to degrees
dtr = pi / 180;              % degrees to radians

% begin main program

clc; home;

fprintf('              program sunsync1 \n\n');

fprintf('   < sun-synchronous orbits - j2 solution > \n\n');

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

   ha = input('please input the mean apogee altitude (kilometers)? ');

   % compute perigee and apogee geocentric radii

   rp = req + hp;
   
   ra = req + ha;

   % compute semimajor axis and eccentricity

   sma = 0.5 * (rp + ra);
   
   ecc = (ra - rp) / (ra + rp);
   
end

% required nodal regression rate (radians/second)

drdt = dtr * (360.0 / 365.2422) / 86400.0;

% Keplerian mean motion (radians/second)

mm = sqrt(mu / (sma * sma * sma));

% orbital semiparameter

slr = sma * (1.0 - ecc * ecc);
   
% inclination initial guess parameter

ipar0 = -(2 / 3) * slr * slr * drdt / (j2 * req * req * mm);

% check for possible solution

if (abs(ipar0) > 1)
    
   clc; home;
   
   disp('no sun-synchronous solution !!');
   
   pause;
   
   inc1 = 0.0;
   
else
    
   % perform iteration

   inc0 = acos(ipar0);

   niter = 0;

   while (1)
       
      niter = niter + 1;

      % perturbed mean motion (radians/second)

      pmm = mm * (1.0 + 1.5 * j2 * req * req * sqrt(1.0 - ecc * ecc) ...
          * (1.0 - 1.5 * sin(inc0) * sin(inc0)) / (slr * slr));

      % updated inclination

      inc1 = acos(-(2 / 3) * slr * slr * drdt / (j2 * req * req * pmm));

      % check for convergence or more than 100 iterations

      if (abs(inc1 - inc0) <= 0.00000001 || niter > 100)
          
         break;
         
      end

      % reset inclination

      inc0 = inc1;
      
   end

end

% print results

clc; home;

fprintf('\n              program sunsync1 \n\n');
   
fprintf('   < sun-synchronous orbits - j2 solution > \n\n');
     
fprintf('mean perigee altitude (kilometers)    %12.4f \n\n', hp);
 
fprintf('mean apogee altitude  (kilometers)    %12.4f \n\n', ha);

fprintf('mean semimajor axis   (kilometers)    %12.4f \n\n', sma);
     
fprintf('mean orbital eccentricity             %12.10f \n\n', ecc);

fprintf('mean orbital inclination (degrees)    %12.4f \n\n', inc1 * rtd);

fprintf('number of iterations                  %12.4f \n\n', niter);


