clc; clear all; close all;

r0=[5662.1 6538.0 3269 ];
v0=[-3.8856 5.1214 -2.2433 ];

h=cross(r0,v0);

i=acos(h(3)/norm(h));

n=cross([0 0 1],h);

Omega=2*pi-acos(n(1)/norm(n));

e=(cross(v0,h)-398600*r0/norm(r0))/398600;

v_r=(r0*v0.')/norm(r0);
e1=((norm(v0)^2-398600/norm(r0))*r0-(r0*v0.')*v0)/398600;
 
omega=acos((n*e.')/(norm(n)*norm(e)));
v_r=r0*v0.'/norm(r0)
theta=acos(e*r0.'/(norm(e)*norm(r0)));
a=norm(h)^2/(398600*(1-norm(e)^2));
