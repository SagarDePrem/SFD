a=140000;
e=.95;
b=a*sqrt(-e^2+1)

%plots the trajectory

x1=linspace(-a,a,1000);
y1= b*sqrt(-x1.^2/a^2+1);
y1a= -b*sqrt(-x1.^2/a^2+1);
plot(x1,y1);
hold on
plot(x1,y1a);
x2=linspace( -a, a,1000);
y2=b*sqrt(-x2.^2/a^2+1);
y2a=-b*sqrt(-x2.^2/a^2+1);
plot(x2,y2)
plot(x2,y2a)

newLim = get(gca,'XLim'); 
newx = linspace(newLim(1), newLim(2),11); 
set(gca,'XTick', newx); 
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
 
 plot(a*e, 0,'.b', 'MarkerSize',20)
% plot(a*e-402000*cos(theta), 402000*sin(theta), '.r', 'MarkerSize',15)
 
 xlabel('X (km)')
ylabel('Y (km)')
title('Trajectory of the meteroid')
