% Works out depths of microcats from pressure data
%
% Required inputs:- 
%      moor - mooring name as string e.g. 'wb_1_200420'
%      pathosnap - path of the osnap directory (ex: ~/Data/osnap)
%
% example:
% ctd_instrdpth2('wb1_1_200420',pathosnap)
%
% Outputs:- outfile -
%           
% Uses the following functions:-
%  rodbload, rodbsave, hms2h

%
%    


function ctd_instrdpth2(moor,pathprocess)
sensor_id = [337];
sensor      = 'microcat';
info_dir    = [pathprocess '/data/moor/proc/', moor, '/'];
info_file   = [moor,'info.dat'];    

infofile   = [info_dir,info_file];

% --- get mooring information from infofile 

[id,sn,z,rcmc1,rcmc2,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:RCMC1:RCMC2:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');


if sensor_id == 337
jj = find(id == 337);
sn_mc = sn(jj);
z_mc = z(jj);
elseif sensor_id==330
    jj = find(id ==330);
    sn_mc = sn(jj)
    z_mc = z(jj);
end
dummy    = -9999

%-----------------------------------------
% --- preprocessing loop -------------------
% ----------------------------------------

for proc = 1 : length(sn_mc)
    if sensor_id==337
  infile  = [info_dir,'microcat/',moor,'_',sprintf('%d',sn_mc(proc)),'.use']
  if exist(infile)  < 1
      infile = [info_dir,'microcat/',moor,'_0',sprintf('%d',sn_mc(proc)),'.use']
  end
    elseif sensor_id == 330
       infile  = [info_dir,'rbr/',moor,'_',sprintf('%d',sn_mc(proc)),'.use']
    end
  if exist(infile)   > 0 

    [idp,sd,st,ed,et,YY,MM,DD,HH,T,C,P] = rodbload(infile,'InstrDepth:Start_Date:Start_Time:End_Date:End_Time:YY:MM:DD:HH:T:C:P');
    
    rodbfile= [infile];
  
 outfile = infile;

 jd       = julian(YY,MM,DD,HH);
 jd0      = jd - jd(1);  
 
Pp=P;
Pp(find(Pp==0))=NaN;
 clf
 plot(jd0,Pp)
hold on
plot(jd0,sw_dpth(Pp,lat),'r')
grid on

mean_pres=nanmean(Pp(10:end-10));
min_p=nanmin(Pp(10:end-10));
depth_mean = sw_dpth(mean_pres,lat);
depth_min = sw_dpth(min_p,lat);

disp(['Mean Pressure:           ', num2str(mean_pres)])
disp(['Minimum Pressure:        ', num2str(min_p)])
disp(['Mean Depth:              ', num2str(depth_mean)])
disp(['Minimum Depth:           ', num2str(depth_min)])
disp(['Instrument header depth: ', num2str(idp)])

replace = input('Do you want to change nominal depth of instrument? ','s');
if strcmp(replace, 'y')
    depth_replace = input('Input new instrument depth: ');
    z_mc(proc) = depth_replace;
    z(jj(proc)) = depth_replace;
    
 %-----------------------------------  
    %--- write data to rodb format -----
    %-----------------------------------

    disp(['writing data to ',outfile]) 
     
    fort = '%4.4d  %3.2d  %3.2d  %9.5f  %8.4f  %8.4f  %4.1f';
    cols = 'YY:MM:DD:HH:T:C:P';
    rodbsave(outfile,...
       'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
       fort,...
       lat,lon,cols,sd,st,sn_mc(proc),mr,wd,z_mc(proc),ed,et,...
       [YY MM DD HH T C P]);
  end
  end
end


if  ~isnan(rcmc1)
       fort2 = '%7d %8d %8d %8d %8d'; 
       cols2 = 'z:instrument:serialnumber:RCMC1:RCMC2';
    rodbsave(infofile,'Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:Columns:WaterDepth:Mooring',...
        fort2,...
        s_t,s_d,e_t,e_d,lat,lon,cols2,wd,mr,...
        [z id sn rcmc1 rcmc2]);
       
else
    
          fort2 = '%7d %8d %8d'; 
       cols2 = 'z:instrument:serialnumber';
    rodbsave(infofile,'Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:Columns:WaterDepth:Mooring',...
        fort2,...
        s_t,s_d,e_t,e_d,lat,lon,cols2,wd,mr,...
        [z id sn]);
       
end