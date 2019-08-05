%% Problem 1
coe_1= coe_from_sv( [6837.432552  1868.795099  1455.480629 ],[-2.294079 6.758849  2.049468 ],398600);
fprintf('\n Angular momentum = %s km^3/s \n', coe_1(1))
fprintf('\n Eccentricity = %s \n', coe_1(2))
fprintf('\n Right ascesion = %s degrees\n', coe_1(3)*180/pi)
fprintf('\n Inclination = %s degrees\n', coe_1(4)*180/pi)
fprintf('\n Argument of Perigee = %s degrees \n', coe_1(5)*180/pi)
fprintf('\n True Anomaly = %s degrees\n', coe_1(6)*180/pi)
fprintf('\n Semi major axis = %s km\n', coe_1(7)) 
%% Problem 2
lst_1 = LST(2016, 12, 17, 5.5, 80.1599639);  % Sriharikota
lst_2 = LST(2007, 12, 21, 10, 144+58/60);    % Melbourne
lst_3 = LST(2006, 02, 15, 3, -43-6/60);       % Rio de Janeiro
fprintf('\n Sidereal time at Sriharikota = %s degrees \n', lst_1) 
fprintf('\n Sidereal time at Melbourne = %s degrees \n', lst_2) 
fprintf('\n Sidereal time at Rio de Janeiro = %s degrees \n', lst_3) 
%% Problem 3
[r,v] = rv_from_observe(988, 4.86, 36, 0.59, 36.6,-0.263, 40,35, 0)
%% Problem 4
r1=[5662.1 6538.0 3269 ];
v1=[-3.8856 5.1214 -2.2433 ];
X=-r1/norm(r1);
Z=-cross(r1,v1);Z=Z/norm(Z);
Y=cross(Z,X);Y=Y/norm(Y);
fprintf('\n Transformation matrix\n') 
T=[X; Y; Z] 