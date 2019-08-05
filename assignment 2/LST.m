function lst = LST(y, m, d, ut, EL)
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
% ������������������������������������������������������������
% Subfunction used in the main body:
% ������������������������������������������������������������
function y = zeroTo360(x)
% �������������������������
%
% This subfunction reduces an angle to the range
% 0 - 360 degrees.
%
% x - The angle (degrees) to be reduced
% y - The reduced value
%
% ------------------------------------------------------------
if (x >= 360)
x = x - fix(x/360)*360;
elseif (x < 0)
x = x - (fix(x/360) - 1)*360;
end
y = x;
return