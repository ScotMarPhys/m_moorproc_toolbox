% function adcp2rodb(moor,varargin)
% Function to convert RDI ADCP data to rodb format.
% raw ADCP data in .mat format as exported by WINADCP
%
% remaps bin depths for different speed of sound as input by user or from
% ctd profile. Does not correct velcoties for speed of sound.

%
% required inputs: moor - mooring name e.g. 'wb2_3_200606'
%                  inpath - if not using standard paths - standarad paths
%                  not setup yet. so required input
% optional inputs: procpath - if not using standard paths for info.dat file
%                  outpath - if not using standard paths
%
% functions called: hms2h.m
%                   julian.m
%                   rodbload.m
%

% Rayner, 10th December 2009 - converted from dvs2rodb_01
% Written for ADCPs with or without pressure sensor
% Houpert Loic, 4/10/16, comment the subfunction for the bin remapping
% using a different sound of speed as it didnt work
%
% Houpert Loic, 19/10/20,  correct bug in the code in case the instrument doesnt have a pressure
% sensor and the depth is fixed. 


function adcp2rodb_01(moor, varargin)

global MOORPROC_G

basedir = MOORPROC_G.moordatadir;
cruise= MOORPROC_G.cruise;

if nargin==0
    help adcp2rodb_01
    return
end

if nargin==1
    pd = moor_inoutpaths('adcp',moor);
else
    pd = varargin{1};
end

operator = MOORPROC_G.operator;

% ----- read infofile / open logfile  ------------------------------------

infovar = 'instrument:serialnumber:z:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth'; 
[id,sn,z,s_t,s_d,e_t,e_d,lat,lon,wd]  =  rodbload(pd.infofile,infovar);

fidlog   = fopen(pd.stage1log,'a');
fprintf(fidlog,'Transformation of ADCP .mat data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));
fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n\n\n',lon);

bg = julian([s_d(:)' hms2h([s_t(:)' 0])]); %start
ed = julian([e_d(:)' hms2h([e_t(:)' 0])]); %end

vec=find((id>=319) & (id <=328)); % Possible ADCP codes - taken from IMP moorings package

serial_nums=sn(vec)

% -------- load data --------------
for i = 1:length(vec)
    fprintf('Processing sn %d',serial_nums(i))
    infile=fullfile(pd.rawpath,sprintf('%d_data.mat',serial_nums(i)));
    
    indep=z(vec(i));
    if length(indep)>1
        indep=indep(1);
    end
    fprintf(fidlog,'infile : %s\n',infile);
    fprintf(fidlog,'ADCP serial number  : %d\n',serial_nums(i));

    all_data=load(infile);
    month=all_data.SerMon; day=all_data.SerDay; year=all_data.SerYear+2000; % year is year since 2000 
    hour=all_data.SerHour; minute=all_data.SerMin; second=all_data.SerSec; hund_seconds=all_data.SerHund;
    pitch=all_data.AnP100thDeg/100; roll=all_data.AnR100thDeg/100; 
    heading=all_data.AnH100thDeg/100;
    % transform heading into degrees North instead of mathematical space
    heading2=90-heading;
    heading2(heading>90)=heading2(heading>90)+360;
    heading=heading2;
    
    t=all_data.AnT100thDeg/100;
    Amp1=all_data.SerEA1cnt; Amp2=all_data.SerEA2cnt; Amp3=all_data.SerEA3cnt; Amp4=all_data.SerEA4cnt;
    Beam1Cor=all_data.SerC1cnt; Beam2Cor=all_data.SerC2cnt; Beam3Cor=all_data.SerC3cnt; Beam4Cor=all_data.SerC4cnt;
    u=all_data.SerEmmpersec/10; v=all_data.SerNmmpersec/10; w=all_data.SerVmmpersec/10;
    err=all_data.SerErmmpersec/10;
    spd=all_data.SerMagmmpersec/10; dir=all_data.SerDir10thDeg/10;
    PG1=all_data.SerPG1; PG2=all_data.SerPG2; PG3=all_data.SerPG3; PG4=all_data.SerPG4;
    
    % create array of latitudes the same size as depth for use in sw_pres
    % and sw_dpth
    % sw_lat=zeros(size(u))+lat; !!!!!!!!!
    sw_lat=lat;
    press_sensor=input('\nDoes the ADCP have a pressure sensor? y/n:  ','s');
    if (strcmp(press_sensor,'y')||strcmp(press_sensor,'Y')||strcmp(press_sensor,'yes')||strcmp(press_sensor,'Yes')||strcmp(press_sensor,'YES'))
        depth=all_data.AnDepthmm/1000;
        % convert depth to pressure using standard seawater routines
        p=sw_pres(depth,sw_lat);
        press_sensor=1;
    else
        correct_depth=eval(['input(''\nThe input depth during instrument setup was ' num2str(all_data.AnDepthmm(1)/1000) 'm\n If you want to use a new depth for bin-mapping and processing,\n enter it now in m (leave empty if not): '')']);
        if ~isempty(correct_depth)
            depth=all_data.AnDepthmm*0+ correct_depth;
        else
            depth=all_data.AnDepthmm/1000;
        end
        p=sw_pres(depth,sw_lat);        
    end
    
    
    dat       = [year month day hms2h(hour,minute,(second+hund_seconds))];
    jd        = julian(dat);
    valI      = find(jd<ed & jd>bg);
    
    % determine if want to only process a subset of bins e.g. if range was
    % limited due to low scatterers then do not need to process many extra
    % bins of bad data
    bins_to_process=input('\nHow many bins out from the sensor head do you want to include?\n Enter number or 0 for all (default = all): ');
    if (isempty(bins_to_process) || bins_to_process==0)
        bins_to_process=length(all_data.SerBins);
    end
        
    % determine bin depths - start and end
    % bin depths as mapped by the ADCP using fixed salinity entered during
    % setup by the user, and the temperature sensor at the transducer
    bin_depths=zeros(2,bins_to_process);
    bin_depths(1,1)=all_data.RDIBin1Mid-all_data.RDIBinSize/2;
    for m=2:bins_to_process
        bin_depths(1,m)=bin_depths(1,m-1)+all_data.RDIBinSize;
    end
    bin_depths(2,:)=bin_depths(1,:)+all_data.RDIBinSize;

    % create matrix of bin depths for saving to rodb file
    bin_depths_mid=bin_depths(1,:)+(bin_depths(2,:)-bin_depths(1,:))/2;
    bin_depths_mid=bin_depths_mid'*ones(1,length(u));
    bin_depths_mid=bin_depths_mid';
    
    % check if ADCP upward or downward facing for use when defining bin depths.
    up_down=input('\n Was the ADCP deployed upward (u) or downward (d) looking? u/d:  ','s');
    if (strcmp(up_down,'u')||strcmp(up_down,'U')||strcmp(up_down,'up')||strcmp(up_down,'UP')||strcmp(up_down,'Up'))
        up_down=-1;
    elseif (strcmp(up_down,'d')||strcmp(up_down,'D')||strcmp(up_down,'down')||strcmp(up_down,'DOWN')||strcmp(up_down,'Down'))
        up_down=1;
    else
        disp('Response not recognised. Stopping.')
        return
    end
    
    remap_bins=input('\nDo you want to remap the bin depths using a different speed of sound? \n   y/n:  ','s');
%     remap_bins ='n';
    if (strcmp(remap_bins,'y')||strcmp(remap_bins,'Y')||strcmp(remap_bins,'yes')||strcmp(remap_bins,'Yes')||strcmp(remap_bins,'YES'))
        remap_bins=1;
        [new_bin_mids, sos_source]=remap_bins_for_sos(bin_depths_mid,p,lat,t); %calls sub-function with series of pressure if ADCP has P sensor
        sw_lat2=zeros(size(new_bin_mids))+lat;
        new_bin_mids=sw_pres(new_bin_mids,sw_lat2); % convert depths back to pressure for saving with rodb files
        new_bin_mids=new_bin_mids';
    end
    

    % ----- save data to rodb -----------------
    
    for j=1:bins_to_process
        
        outfile = fullfile(pd.stage1path,sprintf(pd.stage1form,serial_nums(i),j));

        columns = 'YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:BEAM1SS:BEAM2SS:BEAM3SS:BEAM4SS:BEAM1COR:BEAM2COR:BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP';
        fort =('%4.4d %2.2d %2.2d  %6.4f  %4.2f  %4.2f %4.1f %4.1f %4.1f  %4.2f %4.2f %4.2f  %4.1f  %4.1f  %3.0f %3.0f %3.0f %3.0f   %3.0f %3.0f %3.0f %3.0f  %4.1f  %3.0f %3.0f %3.0f %3.0f');
        if remap_bins==1
            if j==1
                fprintf(fidlog,'Bins remapped for changes in speed of sound using: %s\n',sos_source);
            end
            data = [dat depth+up_down*new_bin_mids(:,j) t u(:,j) v(:,j) w(:,j) heading pitch roll spd(:,j) dir(:,j) Amp1(:,j) Amp2(:,j) Amp3(:,j) Amp4(:,j) Beam1Cor(:,j) Beam2Cor(:,j) Beam3Cor(:,j) Beam4Cor(:,j) err(:,j) PG1(:,j) PG2(:,j) PG3(:,j) PG4(:,j)];
        else
            data = [dat depth+up_down*bin_depths_mid(:,j) t u(:,j) v(:,j) w(:,j) heading pitch roll spd(:,j) dir(:,j) Amp1(:,j) Amp2(:,j) Amp3(:,j) Amp4(:,j) Beam1Cor(:,j) Beam2Cor(:,j) Beam3Cor(:,j) Beam4Cor(:,j) err(:,j) PG1(:,j) PG2(:,j) PG3(:,j) PG4(:,j)];
        end
        infovar = ['Mooring:Start_Time:Start_Date:End_Time:End_Date:Latitude:Longitude:WaterDepth:' ...
                   'Columns:SerialNumber:InstrDepth']; 
        rodbsave(outfile,infovar,fort,moor,s_t,s_d,e_t,e_d,lat,lon,wd,columns,...
                 sn(vec(i)),indep,data);
        eval(['disp(''Data written to ' outfile ''')']);
        fprintf(fidlog,'outfile: %s\n',sprintf(pd.stage1form,serial_nums(i),j));
    
    end
    % -------- generate logfile entries --------------

    sz   =   size(data);

    fprintf(fidlog,'Instrument Target Depth[m]  : %d\n',indep);
    fprintf(fidlog,'Start date and time         : %s \n',datestr(gregorian(jd(1))));
    fprintf(fidlog,'End date and time           : %s \n',datestr(gregorian(jd(end))));
    sampling_rate = round(1./median(diff(jd)));
    ex_samples = round((jd(end)-jd(1))*sampling_rate+1);
    fprintf(fidlog,'Sampling Frequency [per day]: %d \n',sampling_rate);
    fprintf(fidlog,'Number of samples           : %d; expected: %d \n',sz(1),ex_samples);

    for k=1:bins_to_process
        if k==1
            m_hdg = median(heading(valI,:));
            m_pit = median(pitch(valI,:));
            m_rol = median(roll(valI,:));
            
            fprintf(fidlog,'Median heading / pitch / roll [deg]                 : %4.1f  %4.1f  %4.1f\n',m_hdg, m_pit, m_rol);
            if press_sensor==12 %
                m_p = median(p(valI,:));
                fprintf(fidlog,'Median pressure of instrument [dbar]                       : %4.1f\n',m_p);
            end
            m_t = median(t(valI,:));
            fprintf(fidlog,'Median temperature [deg C]              20959         : %4.2f\n',m_t);
        end
        
        m_u = median(u(valI,k));
        m_v = median(v(valI,k));
        m_w = median(w(valI,k));
        m_Beam1ss = median(Amp1(valI,k));
        m_Beam2ss = median(Amp2(valI,k));
        m_Beam3ss = median(Amp3(valI,k));
        m_Beam4ss = median(Amp4(valI,k));
        m_Beam1cor = median(Beam1Cor(valI,k));
        m_Beam2cor = median(Beam2Cor(valI,k));
        m_Beam3cor = median(Beam3Cor(valI,k));
        m_Beam4cor = median(Beam4Cor(valI,k));
        m_Beam1PGP = median(PG1(valI,k));
        m_Beam2PGP = median(PG2(valI,k));
        m_Beam3PGP = median(PG3(valI,k));
        m_Beam4PGP = median(PG4(valI,k));
        m_err = median(err(valI,k));
        m_spd = median(spd(valI,k));
        m_dir = median(dir(valI,k));
        
        
        fprintf(fidlog,'\nBin %d : nominally %3.2f m - %3.2f m from sensor head\n\n', k, bin_depths(1,k), bin_depths(2,k));
        fprintf(fidlog,'Median velocity u / v / w [cm/s]                                 : %4.1f  %4.1f  %4.1f\n',m_u, m_v, m_w);
        fprintf(fidlog,'Median signal strength Beam1 / Beam2 / Beam3 / Beam4 [counts]    : %3.0f  %3.0f  %3.0f  %3.0f\n',m_Beam1ss, m_Beam2ss, m_Beam3ss, m_Beam4ss);
        fprintf(fidlog,'Median correlation Beam1 / Beam2 / Beam3 / Beam4 [counts]        : %3.0f  %3.0f  %3.0f  %3.0f\n',m_Beam1cor, m_Beam2cor, m_Beam3cor, m_Beam4cor);
        fprintf(fidlog,'Median percent good pings Beam1 / Beam2 / Beam3 / Beam4 [counts] : %3.0f  %3.0f  %3.0f  %3.0f\n',m_Beam1PGP, m_Beam2PGP, m_Beam3PGP, m_Beam4PGP);
        fprintf(fidlog,'Median velocity error [cm/s]                                     : %4.1f\n',m_err);
        fprintf(fidlog,'Median speed [cm/s]                                              : %4.1f\n',m_spd);
        fprintf(fidlog,'Median direction [cm/s]                                          : %5.2f\n',m_dir);
    end
end  % loop of serial numbers

display(' ')
display('------------------')
display('Stage 1 completed')
display('------------------')
display(' ')

end
