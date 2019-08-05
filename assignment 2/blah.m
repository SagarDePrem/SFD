clear all;
clc; close all;
Q=[2006 02 15 3 -43-6/60];
y=Q(1);
m=Q(2);
d=Q(3);
ut=Q(4);
EL= Q(5);
digits(10)


j0 = J0(y, m, d);
%...Equation 5.49:
j = (j0 - 2451545)/36525;
%...Equation 5.50:
g0 = 100.4606184 + 36000.77004*j + 0.000387933*j^2- 2.583e-8*j^3;
%...Reduce g0 so it lies in the range 0 - 360 degrees
g0 = zeroTo360(g0);
%...Equation 5.51:
gst = g0 + 360.98564724*ut/24;
%...Equation 5.52:
lst = gst + EL;
%...Reduce lst to the range 0 - 360 degrees:
lst = lst - 360*fix(lst/360);
return