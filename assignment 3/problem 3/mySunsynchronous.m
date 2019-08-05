clear all; clc;
K=3;                    % revisit period
N=45;                   % no of orbits during between two revisits
T=(K/N)*24*3600;        % initial guess for time period
n0=2*pi/T;              % initial guess for mean motion
we = 1.99106e-7;
mu=398600;
a=((T/2/pi)^2*mu)^(1/3);         % initial guess for semi major axis
R=6378.145;
J2=0.00108263;
i=acos(-2/3*(a/R)^2*we/(n0*J2));      % initial guess for inclination

tol=1;
counter=0;
% while tol>0.1
    counter=counter+1
    for counter=1:1
odot=0.75*n0*J2*(R/a)^2*(5*cos(i)^2-1);         % J2 affects argument of perigee
deln=-0.75*n0*J2*(R/a)^2*(3*sin(i)^2-2);        % J2 affects the mean motion
n=deln+n0+odot;                                 % new mean motion
T=(2*pi/n);                                     % new time period
a_new=((T/2/pi)^2*mu)^(1/3);                    % new semi major axis
i=acos(-2/3*(a/R)^2*we/(n0*J2));                % new inclination
tol=abs(a_new-a);
a=a_new;
n0=n;
  end
 h=a-R
 
 ground_track(h,i)          % plot ground track taking altitude and inclination

