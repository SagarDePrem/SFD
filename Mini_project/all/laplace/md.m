
function m = md(y,d) %day is the number of days
if ((y/4)- floor(y/4)==0)
md = [31 29 31 30 31 30 31 31 30 31 30 31]; %only for leap year
else
    md = [31 28 31 30 31 30 31 31 30 31 30 31];
end
m = 1;

while d > md(m)
    
    d = d - md(m);
    m = m+1;
end

m = m;
end