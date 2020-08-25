% function [DD MM YY]=DOY2date(DOY,year)
% function to convery day of year to gregorian format date DD MM YY
% inputs: DOY - day of year e.g. 301
%         year - e.g. 2008
% outputs: matrix of date in gregorian format
% written by Darren Rayner - 15/11/08 - cruise D334

function [DD MM YY]=DOY2date(DOY,year)
if nargin~=2
    disp('Invalid number of inputs to DOY2date function.')
    return
end

days=[0 31 28 31 30 31 30 31 31 30 31 30 31]; %days in each month
if floor(year/4)==year/4 %finds if year is leap year
    disp('Leap Year')
    days(3)=29;
else
    disp('Not Leap Year')
end

YY=year;
cum_days=zeros(size(days));
for i=2:length(days)
    cum_days(i)=days(i)+cum_days(i-1);
end

day1=DOY-cum_days;
a=find(day1>=0);
MM=a(end);
DD=day1(a(end));
