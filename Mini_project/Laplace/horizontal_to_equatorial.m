% Author:           Shiva Iyer
% Date Created:     July 10, 2013
%
% Calculate the right ascension <RA> and declination <dec> for a satellite
% with azimuth <azi>, elevation <ele> for a station with coordinates
% <lat,lon> at UTC time <yy,mm,dd,hh,min,ss>. All angles are in radians.
function [RA,dec] = horizontal_to_equatorial(azi, ele, lat, lon, yy, mm, dd, hh, min, ss)
    ca = cos(azi); sa = sin(azi);
    ce = cos(ele); se = sin(ele);
    cl = cos(lat); sl = sin(lat);

    dec = asin(se*sl + ce*cl*ca);
    cd = cos(dec); sd = sin(dec);

    lmst = mod(GMST(yy, mm, dd, hh, min, ss) + lon, 2*pi);
    ha_x = (se - sd*sl) / (cd*cl);
    ha_y = -sa*ce/cd;
    RA = mod(lmst - atan2(ha_y, ha_x), 2*pi);
end