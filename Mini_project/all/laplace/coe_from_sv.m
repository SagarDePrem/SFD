% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function coe = coe_from_sv(R,V,mu)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 %r0=[6837432.552  1868.795099  1455.480629 ];
%v0=[2.294079 6.758849  2.049468 ];
digits(10);
eps = 1.e-10;
r = norm(R);
v = norm(V);
vr = dot(R,V)/r;
H = cross(R,V);
h = norm(H);
%...Equation 4.7:
incl = acos(H(3)/h);
%...Equation 4.8:
N = cross([0 0 1],H);
n = norm(N);
%...Equation 4.9:
if n ~= 0
RA = acos(N(1)/n);
if N(2) < 0
RA = 2*pi - RA;
end
else
RA = 0;
end
%...Equation 4.10:
E = 1/mu*((v^2 - mu/r)*R - r*vr*V);
e = norm(E);
%...Equation 4.12 (incorporating the case e = 0):
if n ~= 0
if e > eps
w = acos(dot(N,E)/n/e);
if E(3) < 0
w = 2*pi - w;
end
else
w = 0;
end
else
w = 0;
end
%...Equation 4.13a (incorporating the case e = 0):
if e > eps
TA = acos(dot(E,R)/e/r);
if vr < 0
TA = 2*pi - TA;
end
else
cp = cross(N,R);
if cp(3) >= 0
TA = acos(dot(N,R)/n/r);
else
TA = 2*pi - acos(dot(N,R)/n/r);
end
end
%...Equation 4.62 (a < 0 for a hyperbola):
a = h^2/mu/(1 - e^2)
coe = [h e RA incl w TA a];
end %coe_from_sv
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~