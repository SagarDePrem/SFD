clear all; clc;
a0=6378.145+(243+275)*0.5;
r=a0; mu=398600;
k=0.67-(r-6378.145-600)*0.02*70/89;
rho=atmosphere(r-6378.145);
aDot=-k*sqrt(mu*a)*0.015*rho*(2*a/r-1)^1.5;
a1=-k*sqrt(mu*a^2)*0.015*rho/2--k*sqrt(mu*a0^2)*0.015*rho/2