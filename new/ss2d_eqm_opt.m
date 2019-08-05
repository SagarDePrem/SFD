function ydot = ss2d_eqm_opt (~, y)

% two-dimensional solar sail polar equations of motion

% required by ss2d_opt.m

% input

%  t    = non-dimensional simulation time
%  y(1) = radial distance (r)
%  y(2) = radial component of velocity (u)
%  y(3) = tangential component of velocity (v)
%  y(4) = polar angle (radians)

% output

%  ydot(1) = r-dot
%  ydot(2) = u-dot
%  ydot(3) = v-dot
%  ydot(4) = theta-dot

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global b1 b2 b3 acc_srp alpha_wrk

% evaluate equations of motion at current conditions

r = y(1);

u = y(2);

v = y(3);
 
afactor = (acc_srp / r^2) * cos(alpha_wrk);

% r-dot

ydot(1) = u; 

% u-dot

ydot(2) = (v^2 / r) - (1.0 / r^2) + afactor * (b1 + b2 * cos(alpha_wrk)^2 ...
    + b3 * cos(alpha_wrk));

% v-dot

ydot(3) = -(u * v / r) + afactor * sin(alpha_wrk) * (b2 * cos(alpha_wrk) + b3);

% theta-dot

ydot(4) = v / r;


