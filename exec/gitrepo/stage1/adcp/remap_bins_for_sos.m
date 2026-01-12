%% Function to apply CTDcorrection for Speed of sound
% written assuming loading a ctd file from D279 which is in .mat format
% also assumes that no problems with primary temp, conductivity and pressure
% sensors

% outputs remapped bins mid-points in distance from sensor head
function [remapped_bin_mid, sos_source]=remap_bins_for_sos(bin_depths_mid,press_depth,lat,ADCP_t)
    cor_type = input(['To correct the bin depths for speed of sound (sos), would you like to use',...
        '\na) sos derived from ctd profiles',...
        '\nb) fixed sos derived from ADCP temp sensor and fixed salinity given by user (in PSU)',...
        '\nc) fixed sos given by user?   a/b/c: '],'s');
    
    if cor_type=='a'
        load_ctd=input(['\nEnter ctd filename (file in .mat format) with',...
            'full path to use a CTD profile to correct the bin depths\n'],'s');
            
       if exist(load_ctd,'file')
           load(load_ctd);
           sos_source=load_ctd;
           a=who('ctd*');
           if isempty(a)
               disp('CTD file format not recognised. Stopping routine.')
               return
           end
           a=char(a); % to convert from cell
           a=a(4:end); % obtains number of CTD station
           sos=eval(['sw_svel(ctd' a '.salin,ctd' a '.temp*1.00024,ctd' a '.press)']);
           ctd_press=eval(['ctd' a '.press']);
           ctd_salin=eval(['ctd' a '.salin']);
           
           % find SOS at sensor head
           ctd_press2=ctd_press'*ones(1,length(press_depth));
           press_depth2=press_depth*ones(1,length(ctd_press));
           press_depth2=press_depth2';
           %ctd_press2=ctd_press2';
           press_diff=abs(ctd_press2-press_depth2); % bin remapping
            % to remap DY120 ADCP to new depth (instrument were setup with wrong depth)
           for ijk=1:length(press_diff)
            sensor_index(ijk)=find(press_diff(:,ijk) == min(press_diff(:,ijk))); % sensor_index and sos will be a vector when have p sensor
           end
           
           c_m=0;
           while c_m==0
               ctd_measured=input(['\nFor sos at sensor, use sos from ctd profile (c) or as measured by\n'...
                                   'ADCP temp sensor with salinity from ctd profile (m)?   c/m: '],'s');
               if (strcmp(ctd_measured,'c')||strcmp(ctd_measured,'C'))
                   sensor_sos=sos(sensor_index); % sensor_index and sos will be single value when don't have p sensor and using sos from ctd profile
                   c_m=1;
               elseif (strcmp(ctd_measured,'m')||strcmp(ctd_measured,'M'))
                   for ijk=1:length(sensor_index)
                       sensor_salin(ijk)=ctd_salin(sensor_index(ijk));
                   end
                   %sensor_salin=sensor_salin*ones(1,length(ADCP_t));
                   sensor_salin=sensor_salin';
                   %keyboard
                   %press_depth=press_depth*ones(size(sensor_salin)); % make vectors same length
                   adcp_lat=lat*ones(size(sensor_salin)); % make vectors same length
                   press_depth=sw_pres(press_depth,adcp_lat); % convert to pressure from depth
                   sensor_sos=sw_svel(sensor_salin,ADCP_t,press_depth);
                   c_m=2;
               else
                   disp('Enter "c" or "m".')
                   c_m=0;
               end
           end
           % convert ctd_press to ctd_depth in m
           adcp_lat2=lat*ones(size(ctd_press)); % first create vector of same size for latitude
           ctd_depth=sw_dpth(ctd_press,adcp_lat2);
       else
           disp('CTD file does not exist! Stopping routine.')
           return
       end
    elseif cor_type=='b'
        % ask for operator's input of fixed salinity value
         
        use_c=input('What value to use for fixed salinity (in PSU) for speed of sound correction? '); 
        sos_source = ['SOS estimated from fixed salinity of',...
            num2str(use_c),' and ADCP temperature record.'];
        s=use_c*ones(size(ADCP_t)); %creates vector of same length as t filled with use_c value for fixed salinity.
        adcp_lat=lat*ones(size(s)); % make vectors same length
        press_depth=sw_pres(press_depth,adcp_lat); % convert to pressure from depth
        sensor_sos=sw_svel(s,ADCP_t,press_depth);
        
        sos_str = input(['Was the instrument set up with a',...
                    'fixed sos of 1500m/s. Enter 0 to proceed (default) or enter ',...
                    'alternative sos: '],'s');
        if strcmp(sos_str,'0') || isempty(sos_str)
            sos=1500;
        else
            sos=str2num(sos_str);
        end
        
    elseif cor_type=='c'
        sensor_sos=input('What value to use for fixed sos');
        
        sos_source = ['Fixed sos of ',...
            num2str(sensor_sos),'m/s as given by user.'];
        
        sos_str = input(['Was the instrument set up with a',...
                    'fixed sos of 1500m/s. Enter 0 to proceed (default) or enter ',...
                    'alternative sos'],'s');
        if strcmp(sos_str,'0') || isempty(sos_str)
            sos=1500;
        else
            sos=str2num(sos_str);
        end
    end

    
%    remapped_bin_mid=zeros(length(press_depth),length(bin_depths_mid));
%    remapped_bin_mid=remapped_bin_mid*NaN;
%    remapped_bin_start=remapped_bin_mid;
%    remapped_bin_end=remapped_bin_mid;
    bin_depths_mid=bin_depths_mid';
    remapped_bin_mid=NaN*zeros(size(bin_depths_mid));
    remapped_bin_start=remapped_bin_mid;
    remapped_bin_end=remapped_bin_mid;
    remapped_bin_start(1,:)=bin_depths_mid(1,:)-(bin_depths_mid(2,:)-bin_depths_mid(1,:))/2; % remapped bin start depth for bin 1 is same as previously
    
    disp('Remapping bins mid-points in distance from sensor head')
    size_bin_depths_mid=size(bin_depths_mid);
    if length(sensor_sos)==1
        sensor_sos=sensor_sos*ones(1,size_bin_depths_mid(2));
    end

   
    for j=1:size_bin_depths_mid(2) %1:number of samples in adcp timeseries
        for i=1:size_bin_depths_mid(1) %1:number of bins to process
            oldLen=bin_depths_mid(2,1)-bin_depths_mid(1,1);
            if cor_type=='a'

                if i==1
                    a=abs(ctd_depth-bin_depths_mid(i,j));
                    a=find(a==min(a));
                else
                    a=abs(ctd_depth-remapped_bin_mid(i-1,j)+newLen);
                    a=find(a==min(a));
                end
            
                a=a(1); % in case mid-depth was exactly half way between two ctd depths
                newLen=oldLen*(sensor_sos(j)/sos(a));
            else
                newLen=oldLen*(sensor_sos(j)/sos);
            end

            %re-assign bins with new depths
            if i>1
                remapped_bin_start(i,j) = remapped_bin_end(i-1,j);
            end
            remapped_bin_end(i,j)=remapped_bin_start(i,j) + newLen;
            remapped_bin_mid(i,j)=remapped_bin_start(i,j) + newLen/2;
            
        end
    end
end