function [f, g] = ss2d_shoot (x)

% objective function and equality constraints

% simple shooting method

% inputs

%  x(1) = current value of transfer time (objective function)
%  x(2, nsegments) = current values of steering angle alpha

% outputs

%  f = vector of equality constraints and
%      objective function evaluated at x

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global tfactor nsegments alpha_wrk xinitial

% compute duration of each time interval (non-dimensional)

deltat = x(1) / nsegments;

% specify number of differential equations

neq = 4;

% truncation error tolerance

tetol = 1.0e-10;

% initialize initial time

ti = -deltat;

% total non-dimensional time of flight

tof = x(1);

% set initial conditions

yi(1) = xinitial(1);

yi(2) = xinitial(2);

yi(3) = xinitial(3);

yi(4) = xinitial(4);

% step size guess (non-dimensional time)

h = (1200.0 / 86400.0) / tfactor;

% integrate for all segments

for i = 1:1:nsegments
    
    alpha_wrk = x(i + 1);
    
    % increment initial and final times
    
    ti = ti + deltat;
    
    tf = ti + deltat;
    
    % integrate from current ti to tf
    
    yfinal = rkf78('ss2d_eqm_opt', neq, ti, tf, h, tetol, yi);
    
    % reset integration vector
    
    yi = yfinal;
    
    % check for end of simulation
    
    if (tf >= tof)
        
        break;
        
    end
    
end

% objective function (minimize non-dimensional transfer time)

f(1) = x(1);

% compute equality constraints (final state boundary conditions)

f(2) = yfinal(1);

f(3) = yfinal(2);

f(4) = yfinal(3);

% transpose

f = f';

% no derivatives

g = [];
