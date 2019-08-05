%For thise who do not have the Financial Toolbox installed
%Replace the months and date functions with this
function [m,d]=monthcust(date,y)
    if (y/4)- floor(y/4)==0
    monthgaps= [31 29 31 30 31 30 31 31 30 31 30 31]; %only for leap year
    else    
    monthgaps=[31,28,31,30,31,30,31,31,30,31,30,31];
    end
    numdays(1)=31;
    for i=2:length(monthgaps)
        numdays(i)=numdays(i-1)+monthgaps(i)
    end
    i=1
    while(date>numdays(i))
        i=i+1;
    end
    m=i;
    if m>1
    d=date-numdays(m-1)
    else d=date
    end
end