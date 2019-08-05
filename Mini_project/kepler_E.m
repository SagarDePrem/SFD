% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function  kepler_E(e,M)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
 e=0.7079772;
 M=144.225*pi/180;
%...Set an error tolerance:
error = 1.e-8;
%...Select a starting value for E:
if M < pi
E = M + e/2;
else
E = M - e/2;
end
%...Iterate on Equation 3.17 until E is determined to within
%...the error tolerance:
ratio = 1;
counter=1;
while abs(ratio) > error
    fprintf(' iteration = %d\n',counter)
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E))
E = E - ratio;
counter=counter+1;
end
 
%kepler_E
fprintf('\n\n--------------------------------------------------------\n\n')
fprintf('\n Converged eccentric anomaly = %s  \n', E)
fprintf('\n--------------------------------------------------------\n\n')
 