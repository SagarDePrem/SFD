clc;clear all; close all;
r1=[5662.1 6538.0 3269 ];
v1=[-3.8856 5.1214 -2.2433 ];
X=-r1/norm(r1);
Z=-cross(r1,v1);Z=Z/norm(Z);
Y=cross(Z,X);Y=Y/norm(Y);
T=[X; Y; Z]