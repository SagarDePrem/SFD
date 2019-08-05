clc; close all; clear all
a=140000  ;
e=0.95;
m1=5.974e24;
m2=00;
mu=6.6742e-20*(m1+m2)
spec_e=-mu/2/a;
r_p=a*(1-e);
h=sqrt(r_p*mu*(1+e));
 
radVel_theta(e,h); % plots radial velocity against theta
normVel_theta(e,h); % plots normal velocity against theta
pathAngle(e); % plots flight path angle against theta

fprintf('\n\n--------------------------------------------------------\n')
fprintf('\n Specific energy\n')
fprintf(' %s\n',spec_e )
fprintf('\n Specific angular momentum \n')
fprintf(' %s\n', h )
fprintf('\n--------------------------------------------------------\n')
