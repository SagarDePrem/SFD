clear all;
clc; close all;
Q=[2007 12 21 10 144+58/60];
y=Q(1);
m=Q(2);
d=Q(3);
UT=Q(4);
E= Q(5);
digits(10)
J0=367*y-floor(7*(y+floor((m+9)/12))/4)+floor(275*m/9)+d+1721013.5;
  
T0=(J0-2451545)/36525 
theta0=100.4606184+36000.77004*T0+0.000387933*T0^2-2.583*10^-8*T0^3 ;  
theta0=theta0-floor(theta0/360)*360;
thetaG=theta0+360.98564724*UT/24; 
g0=theta0; ut=UT;
theta=thetaG+E;  
