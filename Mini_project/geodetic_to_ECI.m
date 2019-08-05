% Author:           Shiva Iyer
% Date Created:     July 10, 2013
%
% Calculate the ECI vector <ECI_vec> and local mean sidereal time <lmst>
% for a station with coordinates <lat,lon,alt> at time <yy,mm,dd,hh,min,ss>.
% All angles are in radians, and distances are in meters.
function [ECI_vec, lmst] = geodetic_to_ECI(lat, lon, alt, yy, mm, dd, hh, min, ss)
    constants;

    sl = sin(lat);
    e2 = f_E * (2 - f_E);
    p = (SMA_E/sqrt(1 - e2*sl*sl) + alt) * cos(lat);
    q = (SMA_E*(1 - e2)/sqrt(1 - e2*sl*sl) + alt) * sl;

    lmst = mod(GMST(yy, mm, dd, hh, min, ss) + lon, 2*pi);
    ct = cos(lmst); st = sin(lmst);

    ECI_vec = [p*ct; p*st; q];
    ECI_vec(4:6, 1) = [-st; ct; 0]*p*OMEGA_E;
    ECI_vec(7:9, 1) = [-ct; -st; 0]*p*OMEGA_E*OMEGA_E;
end