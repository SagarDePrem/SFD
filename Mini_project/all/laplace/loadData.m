%load data function
%This function is used to load the excel data into matlab
%dataFile is the relative file location of the excel file containing data

function [jd,azh,elevation,y,m,d,UT]=loadData(dataFile)

table=readtable(dataFile);%,'sheet','ang_ri1_BL1_24076');

%% Time 
month = zeros(1,length(table.Year));
jd = zeros(1,length(table.Year));
UT = zeros(1,length(table.Year));

for i=1:length(table.Year)
    month(i) = md(table.Year(i),table.Day(i));
    [jd(i),UT(i)]=JD(table.Year(i),month(i),table.Day(i),table.hh(i),table.mm(i),table.ss(i),table.mss(i));
end
elevation=table.Elevation_deg_Tv;
azh=table.Azimuth_deg_;
y=table.Year;
m=month;
d=table.Day;
    
%SPLICING THE BAD DATA
jd=jd(50:650);
azh=azh(50:650);
elevation=elevation(50:650);
y=y(50:650);
m=m(50:650);
d=d(50:650);
UT=UT(50:650);


    


end

