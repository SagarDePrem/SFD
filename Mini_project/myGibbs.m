clear;clc;close all;
%% inputs
load('p6.mat')    % load respective satellite data 
global  f Re wE
deg    = pi/180;
f      = 1/298.256421867;
Re     = 6378.13655;
wE     = 7.292115e-5;           % angular velocity of earth (rad/s)
mu     = 398600.4418;
eEarth = sqrt(1-(1-f)^2);       % Earth's oblateness factor
phi    = 13.0344722;            % BL2 latitude
H      = 0.83991;               % station height
EL     = 77.5109;                 % BL2 stn longitude
j=size(A);                      % number of readings
for i=1:j-1
   Adot(i)=A(i+1)-A(i);         % Azimuth rate
   adot(i)=a(i+1)-a(i);         % elevation rate
end
for s=1:j
  [m(s),d(s)] = md(Year(s),Day(s));     % month and date
end
ut=hh+mm/60+ss/60/60+mss/1000/60/60;    % calculate UT hours
m=m';d=d';                              % month and date
EL   = 77.5109; %BL2 stn longitude
theta=LST(Year, m, d, ut,EL);           % local sidereal time

%%
e=1; %Just an initialization for the loop
kk=0; %Loop counter
while e>0.001 || isnan(coe(2)) || ierr==1
    %   wanted: e below a certain value, and no NaN values
    %   ierr==1, checks for coplanarity of the vectors chosen. if 1, they are not coplanar
    %  Nsite = Re/sqrt(1-(eEarth*sin(phi))^2); %N in slides. A correction for Geodetic frame
    % phiC = atan(1-eEarth^2*(Nsite/(Nsite+H))*tan(phi));
    kk    = kk+1;  
    pts   = sort(randperm(j(1)-1,3)); %randomly chose 3 points from the dataset
    
    for i=1:3
    theta1(i)  = theta(pts(i));   
    rho1(i)    = rho(pts(i));      %range
    rhoDot1(i) = rhoDot(pts(i));   %range rate
    A1(i)      = A(pts(i));         %Azimuth
    Adot1(i)   = Adot(pts(i));     %Azimuth rate
    a1(i)      = a(pts(i));        %Elevation
    adot1(i)   = adot(pts(i));     %Elevation rate
    end
       
    %Using the function rv_from _observe, the state vectors r and v is
    %being calculated. Three vectors are calculated here
    for i=1:3
        [r(:,i),v(:,i)] = rv_from_observe(rho1(i), rhoDot1(i), A1(i), Adot1(i), a1(i), adot1(i), theta1(i), phi, H);
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
fprintf(' Gibbs Solution:')
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
    fprintf('\n    Period (min)               = %g', T/60)
end
fprintf('\n-----------------------------------------------------\n')