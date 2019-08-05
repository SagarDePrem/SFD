clear all; clc;
%coe = coe_from_sv([5662.1 6538.0 3269 ],[-3.8856 5.1214 -2.2433 ],398600)
%Omega=coe(3); omega=coe(5); i=coe(4); theta=coe(6);
 syms Omega omega theta i
R1=[cos(Omega) sin(Omega) 0; -sin(Omega) cos(Omega) 0; 0 0 1];
R3=[cos(omega+theta)  sin(omega+theta) 0;...
    -sin(omega+theta) cos(omega+theta) 0; 0 0 1];
R2=[1 0 0; 0 cos(i)  sin(i); 0 -sin(i) cos(i)];
fprintf('\n The transformation mattrix is \n')
R=-R3*R2*R1
 
