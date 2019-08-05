% ss2d_opt_64bit.m           July 8, 2014

% two-dimensional solar sail trajectory optimization

% minimize transfer time with piecewise-linear steering

% 64 bit SNOPT algorithm (Mar 17, 2014 version)

% time and state variables

%  t    = simulation time
%  y(1) = radial distance (r)
%  y(2) = radial component of velocity (u)
%  y(3) = tangential component of velocity (v)
%  y(4) = polar angle (theta)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

global rkcoef tfactor acc_srp nsegments

global b1 b2 b3 alpha_wrk xinitial xfinal

% conversion factors - radians to/from degrees

rtd = 180.0d0 / pi;

dtr = pi / 180.0;

% astronomical unit (kilometers)

aunit = 149597870.691d0;

% gravitational constant of the sun (au^3/day^2)

smu = 0.2959122082855912D-03;

% non-dimensional to dimensional time conversion factor

tfactor = sqrt(1.0 / smu);

% initialize rkf78 algorithm

rkcoef = 1;

clc; home;

fprintf('\n\n program ss2d_opt \n');

fprintf('\n Two-dimensional solar sail trajectory analysis \n');

[filename, pathname] = uigetfile('*.dat', 'Please select the input file to read');

[fid, nsegments, achar, b1, b2, b3, iplanet, time_g, time_lb, time_ub, ...
    alpha_g, alpha_lb, alpha_ub] = readss_opt(filename);

% compute srp-to-solar acceleration ratio (au/day^2)

acc_srp = (achar * (86400.0d0^2) / (1000.0d0 * aunit)) / smu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial and final times and states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define state vector at initial time

xinitial(1) = 1.0d0;

xinitial(2) = 0.0d0;

xinitial(3) = 1.0d0;

xinitial(4) = 0.0;

% final conditions

if (iplanet == 1)
    
    % Venus
    
    xfinal(1) = 0.723331;
    
else
    
    % Mars
    
    xfinal(1) = 1.52368d0;
    
end

xfinal(2) = 0.0d0;

xfinal(3) = sqrt(1.0d0 / xfinal(1));

% initial guess for non-dimensional transfer time

xg(1) = time_g / tfactor;

% initial guess for steering angles (radians)

xg(2:nsegments + 1) = alpha_g;

% transpose initial guess

xg = xg';

% upper and lower bounds for non-dimensional transfer time

xlb(1) = time_lb / tfactor;

xub(1) = time_ub / tfactor;

% upper and lower bounds for steering angles (radians)

xlb(2:nsegments + 1) = alpha_lb;

xub(2:nsegments + 1) = alpha_ub;

% transpose bounds

xlb = xlb';

xub = xub';

% define lower and upper bounds on objective function (transfer time)

flow(1) = 0.0d0;

fupp(1) = +Inf;

% define bounds on final state vector equality constraints

flow(2) = xfinal(1);
fupp(2) = xfinal(1);

flow(3) = xfinal(2);
fupp(3) = xfinal(2);

flow(4) = xfinal(3);
fupp(4) = xfinal(3);

flow = flow';

fupp = fupp';

xmul = zeros(26, 1);

xstate = zeros(26, 1);

fmul = zeros(4, 1);

fstate = zeros(4, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve solar sail shooting problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snscreen('on');

[x, f, inform, xmul, fmul] = snopt(xg, xlb, xub, xmul, xstate, ...
    flow, fupp, fmul, fstate, 'ss2d_shoot');

% print results

if (iplanet == 1)
    
    fprintf('\n\n program ss2d_opt - Earth-to-Venus\n');
    
else
    
    fprintf('\n\n program ss2d_opt - Earth-to-Mars\n');
    
end

fprintf ('\n number of segments %2i \n', nsegments);

fprintf ('\n initial state vector \n');

fprintf ('\n  radius                %12.8f (AU)\n', xinitial(1));
fprintf ('\n  radial velocity       %12.8f (AU/day)\n', xinitial(2));
fprintf ('\n  transverse velocity   %12.8f (AU/day)\n', xinitial(3));

fprintf ('\n\n final state vector \n');

fprintf ('\n  radius                %12.8f (AU)\n', f(2));
fprintf ('\n  radial velocity       %12.8f (AU/day)\n', f(3));
fprintf ('\n  transverse velocity   %12.8f (AU/day)\n\n', f(4));

fprintf ('\n  total transfer time   %12.8f (non-dimensional)\n', f(1));
fprintf ('\n  total transfer time   %12.8f days \n\n', tfactor * f(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create control, state and trajectory graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute duration of each time interval

deltat = x(1) / nsegments;

% specify number of differential equations

neq = 4;

% truncation error tolerance

tetol = 1.0e-10;

% initialize initial time

ti = -deltat;

% total time of flight

tof = f(1);

% set initial conditions

yi(1) = xinitial(1);

yi(2) = xinitial(2);

yi(3) = xinitial(3);

yi(4) = xinitial(4);

% load initial graphics data

npts = 1;

x1(npts) = 0.0;

x2(npts) = yi(1);

x3(npts) = yi(2);

x4(npts) = yi(3);

x5(npts) = x(2);

x6(npts) = yi(4);

% step size guess (non-dimensional time)

h = (1200.0 / 86400.0) / tfactor;

% integrate for all segments

for i = 1:1:nsegments
    
    % current steering angle
    
    alpha_wrk = x(i + 1);
    
    % increment initial and final times
    
    ti = ti + deltat;
    
    tf = ti + deltat;
    
    % integrate from current ti to tf
    
    yfinal = rkf78('ss2d_eqm_opt', neq, ti, tf, h, tetol, yi);
    
    % create graphics data
    
    npts = npts + 1;
    
    x1(npts) = tf;
    
    x2(npts) = yfinal(1);
    
    x3(npts) = yfinal(2);
    
    x4(npts) = yfinal(3);
    
    x5(npts) = alpha_wrk;
    
    x6(npts) = yfinal(4);
    
    % reset integration vector
    
    yi = yfinal;
    
    % check for end of simulation
    
    if (tf >= tof)
        
        break;
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% create graphic displays
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);

x(end+1) = x(end);

stairs(tfactor * x1(1:end), rtd * x(2:end));

if (iplanet == 1)
    
    title('Two-dimensional Earth-to-Venus Solar Sail Trajectory', 'FontSize', 16);
    
else
    
    title('Two-dimensional Earth-to-Mars Solar Sail Trajectory', 'FontSize', 16);
    
end

xlabel('mission elapsed time (days)', 'FontSize', 12);

ylabel('steering angle (degrees)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 control.eps

figure(2);

plot(tfactor * x1, x2);

if (iplanet == 1)
    
    title('Two-dimensional Earth-to-Venus Solar Sail Trajectory', 'FontSize', 16);
    
else
    
    title('Two-dimensional Earth-to-Mars Solar Sail Trajectory', 'FontSize', 16);
    
end

xlabel('mission elapsed time (days)', 'FontSize', 12);

ylabel('radial distance (au)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 rdistance.eps

figure(3);

plot(tfactor * x1, x3);

if (iplanet == 1)
    
    title('Two-dimensional Earth-to-Venus Solar Sail Trajectory', 'FontSize', 16);
    
else
    
    title('Two-dimensional Earth-to-Mars Solar Sail Trajectory', 'FontSize', 16);
    
end

xlabel('mission elapsed time (days)', 'FontSize', 12);

ylabel('radial velocity (au/day)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 vradial.eps

figure(4);

plot(tfactor * x1, x4);

if (iplanet == 1)
    
    title('Two-dimensional Earth-to-Venus Solar Sail Trajectory', 'FontSize', 16);
    
else
    
    title('Two-dimensional Earth-to-Mars Solar Sail Trajectory', 'FontSize', 16);
    
end

xlabel('mission elapsed time (days)', 'FontSize', 12);

ylabel('transverse velocity (au/day)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 vtransverse.eps

figure(5);

plot(tfactor * x1, rtd * x6);

if (iplanet == 1)
    
    title('Two-dimensional Earth-to-Venus Solar Sail Trajectory', 'FontSize', 16);
    
else
    
    title('Two-dimensional Earth-to-Mars Solar Sail Trajectory', 'FontSize', 16);
    
end

xlabel('mission elapsed time (days)', 'FontSize', 12);

ylabel('polar angle (degrees)', 'FontSize', 12);

grid;

print -depsc -tiff -r300 polar_angle.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create trajectory display
%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6);

hold on;

if (iplanet == 2)
    
    % Mars is the destination planet
    
    axis([-1.7 1.7 -1.6 1.6]);
    
else
    
    % Venus is the destination planet
    
    axis([-1.2 1.2 -1.1 1.1]);
    
end

plot (0, 0, 'y*');

% create array of polar angles (radians)

t = 0: pi / 50.0: 2.0 * pi;

% plot Earth orbit

plot(xinitial(1) * sin(t), xinitial(1) * cos(t), 'Color', 'b');

% plot destination planet orbit

plot(f(2) * sin(t), f(2) * cos(t), 'Color', 'r');

% plot and label transfer orbit

plot(x2.* cos(x6), x2.*sin(x6), 'Color', 'k');

plot(x2(1) * cos(x6(1)), x2(1) * sin(x6(1)), 'ko');

plot(x2(end) * cos(x6(end)), x2(end) * sin(x6(end)), 'ko');

% label plot and axes

xlabel('X coordinate (AU)', 'FontSize', 12);

ylabel('Y coordinate (AU)', 'FontSize', 12);

if (iplanet == 1)
    
    title('Two-dimensional Earth-to-Venus Solar Sail Trajectory', 'FontSize', 16);
    
else
    
    title('Two-dimensional Earth-to-Mars Solar Sail Trajectory', 'FontSize', 16);
    
end

grid;

axis equal;

zoom on;

print -depsc -tiff -r300 trajectories.eps
