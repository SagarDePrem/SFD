% sunsync3.m        September 2, 2013

% calculates the osculating orbital inclination required
% for a user-defined sun-synchronous orbit

% numerical integration solution

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global dtr rtd suncoef rkcoef 

global req mu smu omega

global raan0 ssdraan oev

global tnode j2 lgrav mgrav tnorbits raan_dot

global jdate0 gst0 isun ccoef scoef

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% astrodynamic and utility constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Earth equatorial radius (kilometers)

req = 6378.14;

% Earth gravitational constant (kilometer^3/second^2)

mu = 398600.4415;

% Sun gravitational constant (kilometer^3/second^2)

smu = 132712440040.9446;

% Earth inertial rotation rate (radians/second)

omega = 7.2921151467e-5;

% angular conversion constants

dtr = pi / 180.0;

rtd = 180.0 / pi;

% required nodal regression rate (radians/second)

ssdraan = dtr * (360.0 / 365.2422) / 86400.0;

% root-finding convergence criterion

rtol = 1.0e-8;

% initialize rkf78

rkcoef = 1;

% initialize sun ephemeris

suncoef = 1;

% read gravity model coefficients

[ccoef, scoef] = readgm('egm96.dat');

j2 = -ccoef(3, 1);

% begin simulation

clc; home;
   
fprintf('                  program sunsync3\n');
   
fprintf('\n< sun-synchronous orbits - integrated solution >\n\n');

% request initial calendar date

[month, day, year] = getdate;

% request initial utc

[utc_hr, utc_min, utc_sec] = get_utc_time;

while(1)
      
   fprintf('\nplease input the simulation duration (nodal periods)\n');
   
   tnorbits = input('? ');
   
   if (tnorbits > 0)
      break;
   end 
   
end

% request degree and order of gravity model

while(1)
    
   fprintf('\nplease input the degree of the Earth gravity model (zonals)\n');
   fprintf('(2 <= zonals <= 18)\n');
   
   lgrav = input('? ');
   
   if (lgrav >= 0)
       
      break;
      
   end
   
end   

while(1) 
    
   fprintf('\nplease input the order of the Earth gravity model (tesserals)\n');
   fprintf('(0 <= tesserals <= 18)\n');
   
   mgrav = input('? ');
   
   if (mgrav >= 0)
       
      break;
      
   end
   
end

while(1)
    
   fprintf('\nwould you like to include the point-mass gravity of the Sun (y = yes, n = no)\n');
   
   yn = lower(input('? ', 's'));
   
   if (yn == 'y' || yn == 'n')
       
      break;
      
   end
   
end

if (yn == 'y')
    
   isun = 1;
   
else
    
   isun = 0;
   
end

% request initial orbital elements

oev = getoe([1;1;1;1;1;0]);

% save initial raan value

raan0 = oev(5);

% and initial true anomaly to nodal crossing

oev(6) = -oev(4);

% compute julian date at 0 hours utc

day = day + utc_hr / 24.0 + utc_min / 1440.0 + utc_sec / 86400.0;

jdate0 = julian(month, day, year);

% compute greenwich sidereal time at 0 hours utc

gst0 = gast4(jdate0, 0.0, 1);

% define a 'bracket' for the orbital inclination

x1 = oev(3) - 5.0 * dtr;

x2 = oev(3) + 5.0 * dtr;

% calculate osculating inclination

clc; home;

fprintf('\n  working ...\n');

[xroot, froot] = brent('ss3func1', x1, x2, rtol);

% print results

clc; home;
   
fprintf('\n                program sunsync3\n');
   
fprintf('\n< sun-synchronous orbits - integrated solution >\n');

fprintf('\nsemimajor axis             %12.6f  kilometers \n', oev(1));

fprintf('\neccentricity               %12.8f  \n', oev(2));

fprintf('\ninclination                %12.6f  degrees \n', xroot * rtd);

fprintf('\nargument of perigee        %12.6f  degrees \n', oev(4) * rtd);

fprintf('\ninitial raan               %12.6f  degrees \n', oev(5) * rtd);

fprintf('\naverage nodal period       %12.6f  minutes \n', tnode / 60.0);

fprintf('\ndesired raan_dot           %12.8f  degrees/day \n', 86400.0 * rtd * ssdraan);

fprintf('\npredicted raan_dot         %12.8f  degrees/day \n', 86400.0 * rtd * raan_dot);

fprintf('\ndegree of gravity model        %2i \n', lgrav);
   
fprintf('\norder of gravity model         %2i \n', mgrav);

fprintf('\nnumber of nodal periods        %3i \n', tnorbits);

if (isun == 1)
    
   fprintf('\nsimulation includes the point-mass gravity of the Sun\n\n');
   
else
    
   fprintf('\nsimulation does not include the point-mass gravity of the Sun\n\n');
   
end
   

   
