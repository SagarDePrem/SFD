% Author:           Shiva Iyer
% Date Created:     July 10, 2013
%
% Calculate the Greenwich Mean Sidereal Time <gmst> at UTC time
% <yy,mm,dd,hh,min,ss>. All angles and GMST are in radians.
function gmst = GMST(yy, mm, dd, hh, min, ss)
    constants;

    c = (julian_date(yy, mm, dd, 0, 0, 0) - JD_2000) / CENTURY;
    t = hh*3600 + min*60 + ss;

    gmst = ((-6.2E-6*c + 0.093104)*c + 8640184.812866)*c + 24110.54841;
    gmst = mod(gmst*(pi/43200) + OMEGA_E*t, 2*pi);
end