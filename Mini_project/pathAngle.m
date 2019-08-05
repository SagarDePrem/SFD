function pathAngle(e)

 

del=pi/360;
n=721;
 
for i=1:n+1
theta(i)=(i-1)*del;
gamma(i)=atan2d(e*sin(theta(i)),(1+e*cos(theta(i))));
end

figure
plot(theta*180/pi,gamma)
 
xlabel(' True anomaly \theta (deg)')
ylabel('Flight path angle \gamma (deg)')
title('Flight path angle vs True anomaly \theta')