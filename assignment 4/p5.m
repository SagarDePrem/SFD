v1=[-5.19425 5.19424 -5.19425];
v2=[-8.58481 1.17067 -2.42304];
beta=0.5*(acos(v1*v2'/(norm(v1)*norm(v2))));
delta=pi-2*beta;
mu=398600;
e=1/cos(beta)
a=mu/norm(v1)^2
rp=(e-1)*mu/norm(v1)^2