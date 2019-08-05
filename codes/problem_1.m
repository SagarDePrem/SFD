clc; close all; clear all
%...Input data:
% Moon:
m1 = 5.974e24  
R = 378.137%1737;
m2 = 1000;
t0 = 0;
tf = 250*60;                         %200*hours;
%...End input data

orbit(m1,m2,R,r0,v0,t0,tf)