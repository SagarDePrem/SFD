clear all;
clc;
r=6378.145;
rb=r+637;
mu=398600;
rc=r+485;
rd=r;
e1=(rb-rc)/(rb+rc);
h1=sqrt(rc*mu*(1+e1*cos(0)));
ra=(h1^2/mu)/(1+e1*cos(150*pi/180));
vna1=h1/ra;
vra1=mu*e1*sin(150*pi/180)/h1;
va1=sqrt(vna1^2+vra1^2);
gamma1=atan(vra1/vna1);

e2=(ra-rd)/(rd-ra*cos(150*pi/180));
h2=sqrt(rd*mu*(1+e2));
vna2=h2/ra;
vra2=mu*e2*sin(150*pi/180)/h2;
va2=sqrt(vna2^2+vra2^2);
gamma2=atan(vra2/vna2);
delGamma=gamma2-gamma1;
delVa=sqrt(va1^2+va2^2-2*va1*va2*cos(delGamma))
phi=atan((vra2-vra1)/(vna2-vna1))*180/pi+180