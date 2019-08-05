%This gives month given a day
function MOD=monthOfDay(year,day)
A=cumsum(eomday( year,1:12));
b=find(A<day);
if(length(b)==0) MOD=1;
else MOD=b(end)+1;
end
