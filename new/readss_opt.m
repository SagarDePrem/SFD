function [fid, nsegments, achar, b1, b2, b3, iplanet, time_g, time_lb, time_ub, ...
    alpha_g, alpha_lb, alpha_ub, jd_guess] = readss_opt(filename)

% read simulation data file

% required by ss2d_opt.m

% input

%  filename = name of simulation data file

% output

%  fid       = file id
%  nsegments = number of trajectory time segments
%  achar     = characteristic acceleration (meters/second^2)
%  b1        = optical force model coefficient b1
%  b2        = optical force model coefficient b2
%  b3        = optical force model coefficient b3
%  iplanet   = target planet (1 = venus, 2 = mars)
%  time_g    = initial guess for mission duration (days)
%  time_lb   = lower bound for mission duration (days)
%  time_ub   = upper bound for mission duration (days)
%  alpha_g   = initial guess for steering angles (radians)
%  alpha_lb  = lower bound for steering angles (radians)
%  alpha_ub  = upper bound for steering angles (radians)
%  jd_guess  = julian date guess for departure

% NOTE: all angular elements are returned in radians

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dtr = pi / 180.0;

% open data file

fid = fopen(filename, 'r');

% check for file open error

if (fid == -1)
    
    clc; home;
    
    fprintf('\n\n error: cannot find this file!!');
    
    pause;
    
    return;
    
end

% read 43 lines of data file

for i = 1:1:43
    
    cline = fgetl(fid);
    
    switch i
        
        case 7
            
            % number of trajectory time segments
            
            nsegments = str2num(cline);
            
        case 10
            
            % characteristic acceleration (meters/second^2)
            
            achar = str2double(cline);
            
        case 13
            
            % optical force model coefficient b1
            
            b1 = str2double(cline);
            
        case 16
            
            % optical force model coefficient b2
            
            b2 = str2double(cline);
            
        case 19
            
            % optical force model coefficient b3
            
            b3 = str2double(cline);
            
        case 22
            
            % target planet
            
            iplanet = str2num(cline);
            
        case 25
            
            % initial guess for mission duration (days)
            
            time_g = str2double(cline);
            
        case 28
            
            % lower bound for mission duration (days)
            
            time_lb = str2double(cline);
            
        case 31
            
            % upper bound for mission duration (days)
            
            time_ub = str2double(cline);
            
        case 34
            
            % initial guess for steering angles
            
            alpha_g = dtr * str2double(cline);
            
        case 37
            
            % lower bound for steering angles
            
            alpha_lb = dtr * str2double(cline);
            
        case 40
            
            % upper bound for steering angles
            
            alpha_ub = dtr * str2double(cline);

        case 43

            % initial guess for departure calendar date

            tl = size(cline);

            ci = strfind(cline, ',');

            % extract month, day and year

            month = str2double(cline(1:ci(1)-1));

            day = str2double(cline(ci(1)+1:ci(2)-1));

            year = str2double(cline(ci(2)+1:tl(2)));

            jd_guess = julian(month, day, year);
    end
    
end

fclose(fid);

