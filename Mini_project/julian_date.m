% Author:           Shiva Iyer
% Date Created:     December 26, 2013
%
% Calculate the Julian Date <jd> corresponding to the
% Gregorian date/time <yy,mm,dd,hh,min,ss>.
function jd = julian_date(yy, mm, dd, hh, min, ss)
	m = fix((mm-14)/12);
	jd = fix(dd-32075+(1461*(yy+4800+m))/4) + ...
		fix(367*(mm-2-12*m)/12) - ...
		fix(3*fix((yy+4900+m)/100)/4);
	jd = jd+(hh-12)/24+min/1440+ss/86400;
end