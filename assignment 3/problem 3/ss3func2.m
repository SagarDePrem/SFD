function fx = ss3func2 (x)

% z component of position vector objective function
    
% required by sunsync3.m

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global tiwrk yi yfinal

% step size guess

h = 60.0;

% truncation error tolerance

tetol = 1.0e-10;

% number of differential equations

neq = 6;

% integrate from tiwrk to requested time = x

yfinal = rkf78('ss3eqm', neq, tiwrk, x, h, tetol, yi);

% current value of z component of position vector

fx = yfinal(3);

