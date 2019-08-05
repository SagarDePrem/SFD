% Author:           Shiva Iyer
% Date Created:     July 10, 2013
%
% Physical and mathematical constants used throughout the code

CENTURY     = 36525; % days
RAD_PER_MIN = (2*pi) / 1440; % rad/min
DEG_TO_RAD  = pi/180; % rad/deg
RAD_TO_DEG  = 180/pi; % deg/rad

JD_1950     = 2433281.5; % JD of 31 December, 1949 0h UTC
JD_2000     = 2451545.0; % JD of 01 January, 2000 12h UTC

SMA_E       = 6378137.0; % m
f_E         = 1/298.257223563; % []
OMEGA_E     = 7.2921150E-5; % rad/s
MU_E        = 3.986004418E14; % m^3/s^2
MU_E_ER     = MU_E/(SMA_E*SMA_E*SMA_E); % (ER)^3/s^2

LIGHT_SPEED = 299792458; % m/s