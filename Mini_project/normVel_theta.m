function normVel_theta(e,h)
%input parameters
G = 6.6742e-20;
m1 = 5.974e24; %moon %5.974e24 - earth;
m2= 1000 ;
 
mu= G*(m1+m2);
del=pi/360;
n=721;
 
for i=1:n+1
theta(i)=(i-1)*del;
v_r(i)=mu*(1+e*cos(theta(i)))/h;
end

figure
plot(theta*180/pi, v_r)
xlabel('True anomaly \theta (deg)')
ylabel('Normal velocity v_n km/s')
title('Normal velocity vs \theta')