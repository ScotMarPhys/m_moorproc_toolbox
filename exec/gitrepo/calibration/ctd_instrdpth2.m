% Works out depths of microcats or aquadops from pressure data, allows selection of a
% modified nominal depth, and writes to file
%
% Required inputs:-
%      moor - mooring name as string e.g. 'wb_1_200420'
%      pathosnap - path of the osnap directory (ex: ~/Data/osnap) that
%      contains data/moor/proc
% Optional input:-
%      sensor_id - 337 for SBE microcat (default) or 330 for RBR or 370 for
%      Nortek
%
% example:
% ctd_instrdpth2('wb1_1_200420',pathosnap)
%
% Outputs:-
%      for microcats, (over)writes rodb format info file for each sensor
%      (if depth changed)
%      For microcats or nortek, overwrites {moor}info.dat for mooring moor
%      (if any sensor depths changed)
%
% Uses the following functions:-
%  rodbload, rodbsave
%  sw_dpth
%


function ctd_instrdpth2(moor,varargin)
global MOORPROC_G

if nargin>2
    sensor_id = varargin{1};
else
    sensor_id = 337;
end
info_dir  = fullfile(MOORPROC_G.moordatadir, 'proc', moor);
infofile  = fullfile(info_dir, [moor 'info.dat']);
if sensor_id == 337
    datadir = fullfile(info_dir, 'microcat');
elseif sensor_id == 330
    datadir = fullfile(info_dir, 'rbr');
elseif sensor_id == 370
    datadir = fullfile(info_dir, 'nor');
end

% --- get mooring information from infofile

[id,sn,z,rcmc1,rcmc2,s_t,s_d,e_t,e_d,lat,lon,wd,mr]  =  rodbload(infofile,'instrument:serialnumber:z:RCMC1:RCMC2:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:Mooring');

% only examine this type of sensor
jj = find(id == sensor_id);
sn_mc = sn(jj);
z_mc = z(jj);

%-----------------------------------------
% --- preprocessing loop -------------------
% ----------------------------------------

changed = 0;
for proc = 1 : length(sn_mc)
    infile  = fullfile(datadir, sprintf('%s_%d.use',moor,sn_mc(proc)));
    if exist(infile, 'file')  < 1
        infile = fullfile(datadir, sprintf('%s_0%d.use',moor,sn_mc(proc)));
    end
    fprintf(1,'infile = %s\n',infile);
    if exist(infile,'file')   > 0
        
        if sensor_id==370
            [idp,sd,st,ed,et,YY,MM,DD,HH,T,P] = rodbload(infile,'InstrDepth:Start_Date:Start_Time:End_Date:End_Time:YY:MM:DD:HH:T:P');
        else
            [idp,sd,st,ed,et,YY,MM,DD,HH,T,C,P] = rodbload(infile,'InstrDepth:Start_Date:Start_Time:End_Date:End_Time:YY:MM:DD:HH:T:C:P');
        end
        
        outfile = infile;
        
        dv = [YY MM DD HH zeros(length(HH),2)];
        dd = datenum(dv); dd = dd-dd(1);
        
        Pp=P;
        Pp(Pp==0)=NaN;
        clf
        plot(dd,Pp,'b',dd,sw_dpth(Pp,lat),'r')
        grid on
        xlabel('days since start')
        ylabel('blue: p (dbar); red z (m)')
        
        %***ylf: cutting off first 9 and last 10 values seems arbitrary, should
        %this be customisable?
        mean_pres = nanmean(Pp(10:end-10));
        min_p = min(Pp(10:end-10));
        depth_mean = sw_dpth(mean_pres,lat);
        depth_min = sw_dpth(min_p,lat);
        
        disp(['Mean Pressure:           ', num2str(mean_pres)])
        disp(['Minimum Pressure:        ', num2str(min_p)])
        disp(['Mean Depth:              ', num2str(depth_mean)])
        disp(['Minimum Depth:           ', num2str(depth_min)])
        disp(['Instrument header depth: ', num2str(idp)])
        
        replace = input('Do you want to change nominal depth of instrument? ','s');
        if strcmp(replace, 'y')
            changed = 1;
            depth_replace = input('Input new instrument depth: ');
            z_mc(proc) = depth_replace;
            z(jj(proc)) = depth_replace;
            
            %-----------------------------------
            %--- write data to rodb format -----
            %-----------------------------------
            
            
            if sensor_id ~= 370
                disp(['writing data to ',outfile])
                fort = '%4.4d  %3.2d  %3.2d  %9.5f  %8.4f  %8.4f  %4.1f';
                cols = 'YY:MM:DD:HH:T:C:P';
                data = [YY MM DD HH T C P];
                rodbsave(outfile,...
                    'Latitude:Longitude:Columns:Start_Date:Start_Time:SerialNumber:Mooring:WaterDepth:Instrdepth:End_Date:End_Time',...
                    fort,...
                    lat,lon,cols,sd,st,sn_mc(proc),mr,wd,z_mc(proc),ed,et,...
                    data);
            end
            
        end
    end
end


if changed
    disp(['writing data to ' infofile])
    %overwrite infofile
    if ~isnan(rcmc1)
        fort2 = '%7d %8d %8d %8d %8d';
        cols2 = 'z:instrument:serialnumber:RCMC1:RCMC2';
        data2 = [z id sn rcmc1 rcmc2];
    else
        fort2 = '%7d %8d %8d';
        cols2 = 'z:instrument:serialnumber';
        data2 = [z id sn];
    end
    rodbsave(infofile,'Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:Columns:WaterDepth:Mooring',...
        fort2,...
        s_t,s_d,e_t,e_d,lat,lon,cols2,wd,mr,data2);
end
