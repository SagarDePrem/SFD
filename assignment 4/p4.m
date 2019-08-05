clear ;clc
Rse=149600000;
Rsv=227.9*10^6;%108200000;
Re=6378.1363;
Rv=6052;
a=0.5*(Rse+Rsv);
mus=132712440018;
vap=sqrt(2*mus*(1/Rse-1/2/a))
ve=sqrt(398600/(Re+400))
vp=(vap^2+2*398600/(400+Re))^0.5
deltaV1=vp-ve;
vv=sqrt(324900/(Rv+1000))
vpe=sqrt(2*mus*(1/Rsv-1/2/a))
vp2=(vpe^2+2*398600/(1000+Rv))^0.5
deltav2=vp2-vv
DELTA_V=  
