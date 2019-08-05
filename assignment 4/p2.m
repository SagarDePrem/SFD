clear all
clc
mu=398600;
a1=15000;
ra2=22000;
rp2=6878;
e2=(ra2-rp2)/(ra2+rp2);
h1=sqrt(398600*15000);
v1=h1/a1;

h2=sqrt(rp2*398600*(1+e2));
theta=acos((h2^2/mu/a1-1)/e2);
v2p=mu*(1+e2*cos(theta))/h2;
v2r=mu*e2*sin(theta)/h2;
v2=(v2p^2+v2r^2)^0.5;
gamma=atan(v2r/v2p);
delV=(v1^2+v2^2-2*v1*v2*cos(gamma))^0.5
alpha=acos((delV^2+v1^2-v2^2)/2/delV/v1)*180/pi
% plot(theta,delV)
% xlabel('True anomaly \theta')
% ylabel('delta-V')
% title('Delta-V required for firing at true anomaly \theta')

% [a b] =min(delV);
% fprintf('The minimum delta-V for the required mission is %s\n',a)
% fprintf('The delta-V must be given at the true anomaly %d degrees of the final orbit\n',theta(b)*180/pi)
% fprintf('The flight path angle is %s\n', gamma(b)*180/pi)
