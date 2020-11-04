function out = save_seacat(outfile,P,T,C,jd,moor,lat,lon,wd,sdate,stime,edate,etime,idepth,serialnumber);
% function out = save_seacat(outfile,P,T,C,jd,moor,lat,lon,wd,sdate,stime,edate,etime,idepth,serialnumber);
%
% save seagauge data to rodb format
%

select = 1; % 1 = P,T // 2 = P,T,C

cols1      = 'YY:MM:DD:HH:P:T';
cols2      = 'YY:MM:DD:HH:P:T:C';

fort1      = '%4.4d  %2.2d  %2.2d  %7.5f   %8.4f %7.4f';
fort2      = '%4.4d  %2.2d  %2.2d  %7.5f   %8.4f %7.4f %7.4f';

time      = gregorian(jd);

data      = [time(:,1:3) hms2h(time(:,4:6)) P(:) T(:)];
if ~isnan(C)
  data = [data C(:)];
  select = 2;
end

switch select

  case 1
   fort = fort1;
   cols = cols1;
  case 2
   fort = fort2;
   cols = cols2;
end  

rodbsave(outfile,...
       'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
         fort,...
         lat,lon,cols,sdate,stime(1:2),serialnumber,moor,wd,idepth,edate,etime(1:2),...
         data);

out =1;
