% function adcp2rodb(moor,varargin)
% Function to convert S55 ADCP data to rodb format.
% raw ADCP data in .mat format as exported by Signature Viewer
%
% This does not work at the moment: remaps bin depths for different speed of sound as input by user or from
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
%
%K Burmeister - converted from RDI to S55 adcp

function adcp2rodb_02(moor, varargin)

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

    in_data=load(infile);
    all_data = in_data.Data;
    all_data = cell2struct(cellfun(@double,struct2cell(all_data),'uni',false),fieldnames(all_data),1);
    [year,month,day,hour,minute,seconds] = datevec(all_data.Average_Time);
    pitch=all_data.Average_Pitch; roll=all_data.Average_Roll; 
    heading=all_data.Average_Heading;

    t=all_data.Average_Temperature;
    Amp1=all_data.Average_AmpBeam1; Amp2=all_data.Average_AmpBeam2; 
    Amp3=all_data.Average_AmpBeam3; Amp4=Amp3*NaN;
    Beam1Cor=all_data.Average_CorBeam1; Beam2Cor=all_data.Average_CorBeam2;
    Beam3Cor=all_data.Average_CorBeam3; Beam4Cor=Beam3Cor*NaN;
    u=all_data.Average_VelEast*100; v=all_data.Average_VelNorth*100; 
    w=all_data.Average_VelUp*100;
    err=all_data.Average_Error;
    spd=err*NaN; dir=err*NaN;
    PG1=err*NaN; PG2=err*NaN; PG3=err*NaN; PG4=err*NaN;
    
    p=all_data.Average_Pressure;
    depth = sw_dpth(p,lat);
    press_sensor=1;
    
    dat       = [year month day hms2h(hour,minute,seconds)];
    jd        = julian(dat);
    valI      = find(jd<ed & jd>bg);
    
    % determine if want to only process a subset of bins e.g. if range was
    % limited due to low scatterers then do not need to process many extra
    % bins of bad data
    Amp1_pro = mean(Amp1);SB1 = find(islocalmin(Amp1_pro),1,'last');
    Amp2_pro = mean(Amp2);SB2 = find(islocalmin(Amp2_pro),1,'last');
    Amp3_pro = mean(Amp3);SB3 = find(islocalmin(Amp3_pro),1,'last');
    SB = min([SB1,SB2,SB3]);
    figure(1),clf
    plot(Amp1_pro),hold on,plot(Amp2_pro),plot(Amp3_pro)
    plot(SB,Amp1_pro(SB),'o')
    ylabel('Temporal mean amplitude [dB]')
    legend('Beam1','Beam2','Beam3')
    title([num2str(SB),' valid bins before surface'])
    bins_to_process=input(['\nAutodetected ', num2str(SB),...
        ' valid bins out from the sensor head.',...
        ' \nDo you want to adjust valid bin number?',...
        ' \nEnter 0 for no (default) or new adjusted number: ']);
    if (isempty(bins_to_process) || bins_to_process==0)
        bins_to_process=SB;
    end
        
    % determine bin depths - start and end
    % bin depths as mapped by the ADCP using fixed salinity entered during
    % setup by the user, and the temperature sensor at the transducer
    bin_depths=zeros(2,bins_to_process);
    blank_dist = in_data.Config.Average_BlankingDistance;
    cell_size = in_data.Config.Average_CellSize/2;
    bin_depths(1,1)=blank_dist+cell_size/2;
    for m=2:bins_to_process
        bin_depths(1,m)=bin_depths(1,m-1)+cell_size;
    end
    bin_depths(2,:)=bin_depths(1,:)+cell_size;

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
    if (strcmp(remap_bins,'y')||strcmp(remap_bins,'Y')||strcmp(remap_bins,'yes')||strcmp(remap_bins,'Yes')||strcmp(remap_bins,'YES'))
        remap_bins=1;
        [new_bin_mids, sos_source]=remap_bins_for_sos(bin_depths_mid,p,lat,t); %calls sub-function with series of pressure if ADCP has P sensor
        sw_lat2=zeros(size(new_bin_mids))+lat;
        new_bin_mids=sw_pres(new_bin_mids,sw_lat2); % convert depths back to pressure for saving with rodb files
        new_bin_mids=new_bin_mids';
    elseif (strcmp(remap_bins,'n')||strcmp(remap_bins,'N')||strcmp(remap_bins,'no')||strcmp(remap_bins,'No')||strcmp(remap_bins,'NO'))
    end
    

    % ----- save data to rodb -----------------
    
    for j=1:bins_to_process
        
        outfile = fullfile(pd.stage1path,sprintf(pd.stage1form,serial_nums(i),j));

        columns = ['YY:MM:DD:HH:Z:T:U:V:W:HDG:PIT:ROL:CS:CD:'...
            'BEAM1SS:BEAM2SS:BEAM3SS:BEAM4SS:BEAM1COR:BEAM2COR:'...
            'BEAM3COR:BEAM4COR:EV:BEAM1PGP:BEAM2PGP:BEAM3PGP:BEAM4PGP'];
        fort =['%4.4d %2.2d %2.2d %6.4f %4.2f %4.2f %4.1f %4.1f %4.1f '...
            '%4.2f %4.2f %4.2f %4.1f %4.1f '...
            '%3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %3.0f %4.1f '...
            '%3.0f %3.0f %3.0f %3.0f'];
        if remap_bins==1
            if j==1
                fprintf(fidlog,'Bins remapped for changes in speed of sound using: %s\n',sos_source);
            end
            data = [dat depth+up_down*new_bin_mids(:,j) t u(:,j) v(:,j) w(:,j) heading pitch roll spd dir Amp1(:,j) Amp2(:,j) Amp3(:,j) Amp4(:,j) Beam1Cor(:,j) Beam2Cor(:,j) Beam3Cor(:,j) Beam4Cor(:,j) err PG1 PG2 PG3 PG4];
        else
            data = [dat depth+up_down*bin_depths_mid(:,j) t u(:,j) v(:,j) w(:,j) heading pitch roll spd dir Amp1(:,j) Amp2(:,j) Amp3(:,j) Amp4(:,j) Beam1Cor(:,j) Beam2Cor(:,j) Beam3Cor(:,j) Beam4Cor(:,j) err PG1 PG2 PG3 PG4];
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
        
        m_u = median(u(valI,k),"omitnan");
        m_v = median(v(valI,k),"omitnan");
        m_w = median(w(valI,k),"omitnan");
        m_Beam1ss = median(Amp1(valI,k),"omitnan");
        m_Beam2ss = median(Amp2(valI,k),"omitnan");
        m_Beam3ss = median(Amp3(valI,k),"omitnan");
        m_Beam4ss = median(Amp4(valI,k),"omitnan");
        m_Beam1cor = median(Beam1Cor(valI,k),"omitnan");
        m_Beam2cor = median(Beam2Cor(valI,k),"omitnan");
        m_Beam3cor = median(Beam3Cor(valI,k),"omitnan");
        m_Beam4cor = median(Beam4Cor(valI,k),"omitnan");

        m_err = median(err(valI),"omitnan");
        m_spd = median(spd(valI),"omitnan");
        m_dir = median(dir(valI),"omitnan");
        
        
        fprintf(fidlog,'\nBin %d : nominally %3.2f m - %3.2f m from sensor head\n\n', k, bin_depths(1,k), bin_depths(2,k));
        fprintf(fidlog,'Median velocity u / v / w [cm/s]                                 : %4.1f  %4.1f  %4.1f\n',m_u, m_v, m_w);
        fprintf(fidlog,'Median signal strength Beam1 / Beam2 / Beam3 / Beam4 [counts]    : %3.0f  %3.0f  %3.0f  %3.0f\n',m_Beam1ss, m_Beam2ss, m_Beam3ss, m_Beam4ss);
        fprintf(fidlog,'Median correlation Beam1 / Beam2 / Beam3 / Beam4 [counts]        : %3.0f  %3.0f  %3.0f  %3.0f\n',m_Beam1cor, m_Beam2cor, m_Beam3cor, m_Beam4cor);
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




