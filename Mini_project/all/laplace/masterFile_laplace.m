%This script file calls other functions and does Orbit determination using
%laplace method.
%Referred paper:
%http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1991SvA....35..428K&defaultprint=YES&page_ind=0&filetype=.pdf
clear all;
close all;

% Re - equatorial radius of the earth (km)
% f - earth's flattening factor
% wE - angular velocity of the earth (rad/s)
global f Re wE deg
deg = pi/180;
Re =  6378.137;
f = 0.003353;
wE = 7.2921159* 10^-5 ;
omega = [0 0 wE];
muEarth=398600;

%% DataDetails
H=839.91;%m Bangalore elevation above sea
TrackingStationNo=1;
longi=360-282.48831;
lati=13.0344722; %geodetic(normal)
dataFile='as1_alldata_laplace.xls';

%%
[jd,azh,elevation,y,m,d,UT]=loadData(dataFile);
lst=zeros(1,length(jd)); %local sidereal time
for i=1:length(jd)
    lst(i)=LST(y(i), m(i), d(i), UT(i), longi);
    [RA(i),dec(i)]=topoToEci(azh(i),elevation(i),lst(i),lati,H/1000);
end
% Cutting the curve to make the samples better
% Only JD, elevation and azhimuth are important from now on
RA=RA';
dec=dec';
fittedRA = spline(UT'*3600,RA) ;%Function of sampke no; not the jd values
RAdot = fnder(fittedRA,1);
RAdot = ppval(RAdot,UT*3600);
RAdot = smooth(RAdot,0.25,'loess');
%plot(ppval(fittedRA,UT*3600))

fitteddec= spline(UT'*3600,dec); %Function of sampke no; not the jd values
decdot = fnder(fitteddec,1);
decdot = ppval(decdot,UT*3600);
decdot = smooth(decdot,0.25,'loess');
%figure
%plot(ppval(fitteddec,UT*3600))
%plot(dec)

%All the future calculations are on new fit values of jd, elevation and
%azhimut

%% Differentiation polynomial
% Dequation=[cos(fittedRA*pi/180).*cos(fitteddec*pi/180);
%     sin(fittedRA*pi/180).*cos(fitteddec*pi/180);
%     sin(fitteddec*pi/180)];% Input is sample number
% DdotEquation=[polyder(Dequation(1,:));
%                 polyder(Dequation(2,:));
%                 polyder(Dequation(3,:))]/jd_samples*1/(24*60*60);
% DdotdotEquation=[polyder(DdotEquation(1,:));
%                 polyder(DdotEquation(2,:));
%                 polyder(DdotEquation(3,:))]/jd_samples*1/(24*60*60);
%

%%
%D vector calculations/substitutions
for i = 1:length(UT)
    D(:,i)=[cosd(RA(i)).*cosd(dec(i));
        sind(RA(i)).*cosd(dec(i));
        sind(dec(i))];
    
    Ddot(:,i) = [-cosd(RA(i)).*sind(dec(i)).*decdot(i) - sind(RA(i)).*cosd(dec(i)).*RAdot(i);
       cosd(RA(i)).*cosd(dec(i)).*RAdot(i) - sind(dec(i)).*sind(RA(i)).*decdot(i);
        cosd(dec(i)).*decdot(i)]*deg;
end
fittedD1 = spline(UT'*3600,Ddot(1,:)');
fittedD2 = spline(UT'*3600,Ddot(2,:)');
fittedD3 = spline(UT'*3600,Ddot(3,:)');
% Ddot1 = ppval(fnder(fittedD1,1),UT*3600);
% Ddot2 = ppval(fnder(fittedD2,1),UT*3600);
% Ddot3 = ppval(fnder(fittedD3,1),UT*3600);
% Ddot1 = smooth(Ddot1,0.25,'loess')';
% Ddot2 = smooth(Ddot2,0.25,'loess')';
% Ddot3 = smooth(Ddot3,0.25,'loess')';
% Ddot = [Ddot1;Ddot2;Ddot3];
Ddotdot1 = ppval(fnder(fittedD1,1),UT*3600);
Ddotdot2 = ppval(fnder(fittedD2,1),UT*3600);
Ddotdot3 = ppval(fnder(fittedD3,1),UT*3600);
Ddotdot1 = smooth(Ddotdot1,0.3,'lowess')';
Ddotdot2 = smooth(Ddotdot2,0.3,'lowess')';
Ddotdot3 = smooth(Ddotdot3,0.3,'lowess')';
Ddotdot = [Ddotdot1;Ddotdot2;Ddotdot3];
%Ddotdot = diff(Ddot,1,2);


%% R calculation

for i=1:length(jd)
    R(:,i)=Rvector(lst(i),lati,H/1000);
    Rdot(:,i) = cross(omega,R(:,i));
    Rdotdot(:,i) = cross(omega,Rdot(:,i));
end

%% Finding P and Q
%r=R+rho D
% rho=P-Q/r^3
coe = 0.02*ones(7,1);
coemin = 10;
amax = 0;
for sampChos=1:length(jd)
    Rs=R(:,sampChos);% R sample of choice
    Ds=D(:,sampChos);% R sample of choice
    Rsdot=Rdot(:,sampChos);
    Rsdotdot=Rdotdot(:,sampChos);
    Dsdot=Ddot(:,sampChos);
    Dsdotdot= Ddotdot(:,sampChos);
    
    D1 = det([Ds,Dsdot,Rsdotdot]);
    D2 = det([Ds,Dsdot,Rs]);
    Det =2*det([Ds,Dsdot,Dsdotdot]);
    F = -2*D1/Det;
    Q =-2*muEarth*D2/Det;
    %     P=-dot(Rsdotdot,cross(Dsdot,Ds)  )/(dot(Dsdotdot,cross(Dsdot,Ds)));
    %     Q=muEarth*( dot(Rs,cross(Dsdot,Ds))   )/(dot(Dsdotdot,cross(Dsdot,Ds)));
         Pdash=-dot(  Rs, cross(Rsdotdot,Ds  )  )...
             /2/dot( Rs, cross(Dsdot,Ds )   );
         Qdash=dot(  Rs, cross(Dsdotdot,Ds  )  )...
             /2/dot( Rs, cross(Dsdot,Ds )   );
    
    %% Have to solve and find rho numerically now using r interms of rho equation
    %     syms rho;
    %     syms func;
    %     func =  Q^2-(rho-P)^2*((norm(Rs))^2 + rho^2 + 2*rho*dot(Ds,Rs))^3 ;
    %     func = vpa(expand(func));
    %     cffs = double(coeffs(func,rho));
    %     %assume(rho, 'real')
    %     Rho=max(roots(cffs));
    
    r = roots(fliplr([Q^2 0 0 2*F*Q+2*Q*dot(Ds,Rs) 0 0 F^2+2*F*dot(Ds,Rs)+norm(Rs)^2 0 -1]));
    r = max(r(imag(r)==0));
    Rho = -2*D1/Det - 2*muEarth*D2/(r^3*Det);
    %% Finding Pos and Velocity in ECI
    pos = Rs+Rho*Ds;
    vel = Rsdot+(Pdash-Qdash*Rho)*Ds+Rho*Dsdot;
    
    %% Keplerian elements finding
    coe = coe_from_sv(pos,vel,muEarth);
    if(coe(2)<coemin)
        coemin = coe(2);
        amax = coe(7);
    end
end
coemin






