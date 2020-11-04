% function [D M Y hh mm ss]=sec2DMYhhmmss(time,year)
% function to convert seconds from start of year to separate values for
%       D,M,Y,hh,mm and ss from seconds from start of year given year
% Required inputs = time (seconds from start of year as with pstar format files)
%                   year (year at start of record - not stored with time variable                  variable)
% written by Darren Rayner 13/11/08 cruise D334

function [D M Y hh mm ss]=sec2DMYhhmmss(time,year)
if nargin~=2
    disp('Invalid number of inputs to sec2DMYhhmmss function.')
    return
end

secs_from_start=time;
mins_from_start=time/60;
hours_from_start=time/3600;
days_from_start=time/(3600*24);

floor_days=floor(days_from_start)+1; % 1st day of year is 1 not 0 so add 1.
days=[0 31 28 31 30 31 30 31 31 30 31 30 31];
if floor(year/4)==year/4 %finds if year is leap year
    disp('Leap Year')
    days(3)=29;
else
    disp('Not Leap Year')
end
D=zeros(size(time)); M=D; Y=D; hh=D; mm=D; ss=D;
Y=Y+year; %currently no provision for files that pass from December in January

for i=1:12
    a=find(and((floor_days<sum(days(1:i+1))), (floor_days>sum(days(1:i)))));
    D(a)=floor_days(a)-sum(days(1:i));
    M(a)=i;
    hh(a)=floor(hours_from_start(a)-(floor(days_from_start(a))*24));
    mm(a)=floor(mins_from_start(a)-(floor(hours_from_start(a))*60));
    ss(a)=floor(secs_from_start(a)-(floor(mins_from_start(a))*60));
end
