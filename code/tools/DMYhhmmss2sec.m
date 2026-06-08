% function secs=DMYhhmmss2sec(gregorian_date);
% function to convert a gregorian date to number of seconds from start of
% year
% Inputs: gregorian_date as single vector [DD MM YYYY hh mm ss]
%   e.g. [13 11 2008 10 52 00]
%
% written by Darren Rayner 13/11/08 - cruise D334

function secs=DMYhhmmss2sec(gregorian_date)

days=[0 31 28 31 30 31 30 31 31 30 31 30 31];
if floor(gregorian_date(3)/4)==gregorian_date(3)/4 %finds if year is leap year
    disp('Leap Year')
    days(3)=29;
else
    disp('Not Leap Year')
end

days_since_start=sum(days(1:(gregorian_date(2)))) + (gregorian_date(1)-1);
hours_since_start=days_since_start*24+gregorian_date(4);
minutes_since_start=hours_since_start*60+gregorian_date(5);
secs=minutes_since_start*60+gregorian_date(6);