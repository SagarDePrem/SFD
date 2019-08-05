function [jd_depart, tsynodic] = departure(jd_guess, transfer_time, transfer_angle)

% estimate departure calendar date

% input

%  jd_guess       = departure julian date initial guess
%  transfer_time  = Earth-to-Mars transfer time (days)
%  transfer_angle = Earth-to-Mars heliocentric transfer angle (degrees)

% output

%  jd_depart = departure julian date
%  tsynodic  = Earth/Mars synodic period (days)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global aunit xmu

pi2 = 2.0 * pi;

% time arguments

t = (jd_guess - 2451545.0) / 36525.0;

t2 = t * t;

% celestial longitude of the Earth and Mars at jd_guess (degrees)

theta_earth =  mod(100.466449 + 35999.3728519 * t - 0.00000568 * t2, 360.0);

theta_mars =  mod(355.433275 + 19140.2993313 * t + 0.00000261 * t2, 360.0);

% orbital period of the Earth and Mars (days)

tau_earth = pi2 * sqrt(aunit^3 / xmu) / 86400.0;

tau_mars = pi2 * sqrt((aunit * 1.52368)^3 / xmu) / 86400.0;

% orbital rate of the Earth and Mars (degrees/day)

omega_earth = 360.0 / tau_earth;

omega_mars = 360.0 / tau_mars;

% compute change to initial guess (days)

delta_t = (theta_mars + omega_mars * transfer_time - theta_earth - ...
    transfer_angle) / (omega_earth - omega_mars);

% departure julian date

jd_depart = jd_guess + abs(delta_t);

% synodic period (days)

tsynodic = 360.0 / abs(omega_earth - omega_mars);
