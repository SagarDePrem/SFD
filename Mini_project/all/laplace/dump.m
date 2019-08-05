%This file is solely used for testing purposes. it is highly gen
function MOD=monthOfDay(year,day)
A=cumsum(eomday( year,1:12))
b=find(A<day)
if(length(b)==0) MOD=1
else MOD=b(end)+1
end
