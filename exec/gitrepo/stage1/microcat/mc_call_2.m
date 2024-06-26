% mc_call_2 is a script that performs stage1 processing
% on microcat data.  It converts microcat data from raw to rodb
% format for an entire mooring.
%
% It calls microcat2rodb (to convert microcat_data), rodbload.m,
% timeaxis.m, auto_filt.m, julian.m 

% 27/10/09 DR added functionality for .cnv files
% 21/03/10 ZBS made all paths rely on a basedir variable,
%          added a descriptive header, and modified for oceanus 459
% 08/11/12 DR changed write mode of log file so appends instead of
%          overwriting
% 12/11/12 DR added check on whether infile exists and clearly highlights
%          it for the operator with the option of stopping the routine.
% This will need functionality for oxygen MicroCATs added prior to the
% recovery of instruments on Autumn 2015 cruise. microcat2rodb_3.m already
% has the functionality and mc_call_caldip_jc103 too so use them as a
% template.

close all
global MOORPROC_G
clearvars -except MOORPROC_G 

% only mooring name needs to be modified, rest set in MOORPROC_G by
% startup{cruise}.m

moor = input('mooring deployment (e.g. ebh2_15_2022) to process:   ','s');
%moor            = 'ebh2_15_2022';
ii = strfind(moor,'_');
YEAR = str2double(moor(ii(end)+1:ii(end)+4)); % year of the first measurement
dateoffset = 0;

cruise          = MOORPROC_G.cruise;
operator        = MOORPROC_G.operator;

%get file paths
pd = moor_inoutpaths('microcat',moor);
out_ext  = '.raw';

% --- get mooring information from infofile ---

[id,sn,lat,lon] = rodbload(pd.infofile,'id:sn:latitude:longitude');
if isempty(id) || isnan(id)
  [id,sn,lat,lon] = ...
      rodbload(pd.infofile,'instrument:serialnumber:latitude:longitude');
end

% --- vector of serial numbers ---
ii = find(id >= 335 & id <=337);
vec = sn(ii);

% --- write header info to log file ---
if ~exist(pd.stage1path,'dir')
    mkdir(pd.stage1path)
end
fidlog = fopen(pd.stage1log,'a');
fprintf(fidlog,'Transformation of ascii data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));

fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n',lon);


% --- loop through each instrument on the mooring ---
for i = 1:length(vec)
    fprintf(fidlog,'\n\n');
    disp(['processing microcat serial number ' num2str(vec(i)) ])

    infiles = {sprintf('%4.4d%s',vec(i),'rec2.asc')
        sprintf('%4.4d%s',vec(i),'rec.asc')
        sprintf('%3.3d%s',vec(i),'rec.asc')
        sprintf('%4.4d%s',vec(i),'REC.asc')
        sprintf('%3.3d%s',vec(i),'REC.asc')
        fullfile('data',sprintf('%4.4d%s',vec(i),'_2.asc'))
        fullfile('data',sprintf('%4.4d%s',vec(i),'.asc'))
        sprintf('%4.4d%s',vec(i),'.asc')
        sprintf('%4.4d%s',vec(i),'_data.cnv')
        sprintf('%4.4d%s',vec(i),'_Data.cnv')
        sprintf('%4.4d%s',vec(i),'_data.asc')
        };
    for n = 1:length(infiles)
        infile = fullfile(pd.rawpath,infiles{n});
        if exist(infile,'file')
            datfileinfo = dir(infile);
            if datfileinfo.bytes>0
               %found it
               break
            end
        end
    end
    infile = fullfile(pd.rawpath,infiles{n});
    if ~exist(infile,'file')
        ssname = regexp(moor,'_','split');
        mcatfname = [ssname{1} '_' sprintf('%4.4d',vec(i)) '_' ...
            num2str(YEAR+1,'%4.0f') '.cnv'];
        infile = fullfile(pd.rawpath,mcatfname);
    end
    stopped=0;
    % if still doesn't exist, flag to operator
    if ~exist(infile,'file')
        disp(' ')
        disp(['WARNING: No data file found for serial number ' num2str(vec(i))])
        want_to_continue=input('Do you wish to continue processing the rest of the mooring? (y/n):- ','s');
        if strncmpi(want_to_continue,'y',1)==0
            disp('Stopping!')
            stopped=1;
            fprintf(fidlog,['\n Routine stopped due to missing datafile for serial number: ' num2str(vec(i))]);
            break
        end
    else
        outfile = fullfile(pd.stage1path,sprintf(pd.stage1form,vec(i)));
        microcat2rodb(infile,outfile,pd.infofile,fidlog,'y',dateoffset);
        disp('Press any key to continue: ')
        pause
    end
end

if stopped==0
    comment = input('Enter additional comment to be save in Log file: ','s'); 
    if ~isempty(comment)
        fprintf(fidlog,'\n COMMENT:\n %s',comment);
    end
end
fclose(fidlog);
