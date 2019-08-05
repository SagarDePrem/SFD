clear;
clc;
close all;
global  f Re wE

deg    = pi/180;
f      = 1/298.256421867;
Re     = 6378.13655;
wE     = 7.292115e-5; %angular velocity of earth (rad/s)
mu     = 398600.4418;
eEarth = sqrt(1-(1-f)^2); %Earth's oblateness factor

myData = xlsread('myDataFile.xlsx', 'all'); %Reading data from excel file
%%
%From the data, change the second parameter of the array.
% e.g, Days in the data is in third column. Hence, myData(1,3) ..Using 1 to
% just take the first value. Similarly for other data
Days = myData(1,3); %Day of year
Y = myData(1,2); %Year
 
[M,D]=md(Y,Days)
%%
e=1; %Just an initialization for the loop
kk=0; %Loop counter
while e>0.001 || isnan(coe(2)) || ierr==1
    %   wanted: e below a certain value, and no NaN values
    %   ierr==1, checks for coplanarity of the vectors chosen. if 1, they are not coplanar
    kk    = kk+1
    pts   = sort(randperm(size(myData,1),3)); %randomly chose 3 points from the dataset
    UT    = myData(pts,4)+myData(pts,5)/60+myData(pts,6)/3600+myData(pts,7)/3600000; %calculate UT in hrs
    lon   = 77.5109; %BL2 stn longitude
    theta = LST(Y, M, D, UT, lon); %function to calc Local Sidereal Time
    
    phi   = deg2rad(13.0344722);  %BL2 latitude
    theta = deg2rad(theta);
    H     = 0.83991; %station height
    Nsite = Re/sqrt(1-(eEarth*sin(phi))^2); %N in slides. A correction for Geodetic frame
    %phiC = atan(1-eEarth^2*(Nsite/(Nsite+H))*tan(phi));
    
    %Columns from the data are assigned
    rho    = myData(pts,9);  %range
    rhodot = myData(pts,30); %range rate
    A      = myData(pts,28); %Azimuth
    Adot   = myData(pts,31); %Azimuth rate
    a      = myData(pts,29); %Elevation
    adot   = myData(pts,32); %Elevation rate
    
    % %     For MT2 sat
    % rho    = myData(pts,9);
    % rhodot = myData(pts,23);
    % A      = myData(pts,21);
    % Adot   = myData(pts,24);
    % a      = myData(pts,22);
    % adot   = myData(pts,25);
    
    %Using the function rv_from _observe, the state vectors r and v is
    %being calculated. Three vectors are calculated here
    for i=1:3
        [r(:,i),v(:,i)] = rv_from_observe(rho(i), rhodot(i), A(i), Adot(i), a(i), adot(i), rad2deg(theta(i)), rad2deg(phi), H);
    end
    %%
    %Gibbs
    %function returns the velocity vector for the central vector and ierr
    %returns about the coplanarity of the three vectors. ierr=1 implies non
    %coplanar vectors and another set of vectors must be chosen
    [v2, ierr] = gibbs(r(:,1), r(:,2), r(:,3));
    
    %...If the vectors r1, r2, r3, are not coplanar, abort:
    if ierr == 1
        fprintf('\n  These vectors are not coplanar.\n\n')
        %return   %Uncomment return to end loop when non coplanar vectors
        %are encountered
    end
    
    %coe_from_sv returns the orbital coefficients from the state vectors
    coe  = coe_from_sv(r(:,2),v2,mu);
    
    %eccentricity from the returned coe
    e    = coe(2);
    
end

%%
%PRINTING DATA
h    = coe(1);
RA   = coe(3);
incl = coe(4);
w    = coe(5);
TA   = coe(6);
a    = coe(7);

%...Output the results to the command window:
fprintf(' Solution:')
fprintf('\n');
fprintf('\n  v2 (km/s) = [%g  %g  %g]', v2(1), v2(2), v2(3))
fprintf('\n\n  Orbital elements:');
fprintf('\n    Angular momentum (km^2/s)  = %g', h)
fprintf('\n    Eccentricity               = %g', e)
fprintf('\n    Inclination (deg)          = %g', incl/deg)
fprintf('\n    RA of ascending node (deg) = %g', RA/deg)
fprintf('\n    Argument of perigee (deg)  = %g', w/deg)
fprintf('\n    True anomaly (deg)         = %g', TA/deg)
fprintf('\n    Semimajor axis (km)        = %g', a)
%...If the orbit is an ellipse, output the period:
if e < 1
    T = 2*pi/sqrt(mu)*coe(7)^1.5;
    fprintf('\n    Period (s)                 = %g', T)
end
fprintf('\n-----------------------------------------------------\n')