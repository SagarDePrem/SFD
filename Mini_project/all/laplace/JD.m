%This fuctnion is used to find julian day(An absolute timing system wrt greenwich)
function [jd,UT] = JD(year, month, day,hh,mm,ss,mss)
j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) ...
+ fix(275*month/9) + day + 1721013.5;
UT=hh+mm/60+ss/60/60+mss/1000/60/60;
jd=j0+UT/24;
end

