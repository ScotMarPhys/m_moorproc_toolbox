% mc_call_2_dy053 is a script that performs stage1 processing
% on microcat data.  It converts microcat data from raw to rodb
% format for an entire mooring.
%
% It calls microcat2rodb_4 (to convert microcat_data), rodbload.m,
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
clearvars  -except pathosnap 

% -----------------------------------------------------------------
% --- This is the information that needs to be modified for -------
% --- different users, directory trees, and moorings --------------
% -----------------------------------------------------------------

cruise = 'msm74';
operator = 'nph';
% moor = 'nocm1_03_2016';
% moor = 'nocm2_03_2016';
% moor = 'nocm3_03_2016';
moor = 'nocm4_03_2016';
% moor = 'nocm5_03_2016';

dateoffset = 2016; % year of the first measurement

if exist('pathosnap','var')
    basedir = [pathosnap '/'];
else
%     basedir = '/home/mstar/osnap/';
    basedir = '/Users/ukosnap/Documents/aaaaMSM74/new_mooring_processing/moorings/';
end


%------------------------------------------

% --- set paths for data input and output ---

addpath(genpath([basedir 'data/exec/' cruise '/']));
%addpath(genpath([basedir 'moor']));

% inpath   = [basedir 'data/moor/raw/' cruise '/microcat/' inputdir];
inpath   = [basedir 'data/moor/raw/' cruise '/microcat/'];
outpath  = [basedir 'data/moor/proc/' moor '/microcat/'];
infofile = [basedir 'data/moor/proc/' moor '/' moor 'info.dat'];

out_ext  = '.raw';

% --- get mooring information from infofile ---

[id,sn,lat,lon] = rodbload(infofile,'id:sn:latitude:longitude');
if isempty(id) | isnan(id)
  [id,sn,lat,lon] = ...
      rodbload(infofile,'instrument:serialnumber:latitude:longitude');
end

% --- vector of serial numbers ---
ii = find(id >= 333 & id <=337);
vec = sn(ii);

% --- write header info to log file ---
fidlog = fopen([outpath,'stage1_log'],'a');
fprintf(fidlog,'Transformation of ascii data to rodb format \n');
fprintf(fidlog,'Processing carried out by %s at %s\n\n\n',operator,datestr(clock));

fprintf(fidlog,'Mooring   %s \n',moor);
fprintf(fidlog,'Latitude  %6.3f \n',lat);
fprintf(fidlog,'Longitude %6.3f \n',lon);


% --- loop through each instrument on the mooring ---
for i = 1:length(vec),
  fprintf(fidlog,'\n\n');
  disp(['processing microcat serial number ' num2str(vec(i)) ])
  
  infile = [inpath,sprintf('%4.4d',vec(i)),'rec2.asc'];

  if exist(infile) ~= 2
    infile = [inpath,sprintf('%4.4d',vec(i)),'rec.asc'];
  end
  if exist(infile) ~= 2
    infile = [inpath,sprintf('%3.3d',vec(i)),'rec.asc'];
  end
  if exist(infile) ~= 2
    infile = [inpath,sprintf('%4.4d',vec(i)),'REC.asc'];
  end
  if exist(infile) ~= 2
    infile = [inpath,sprintf('%3.3d',vec(i)),'REC.asc'];
  end

  if exist(infile) ~= 2
    infile = [inpath,'data',sprintf('%4.4d',vec(i)),'_2.asc'];
  end

  if exist(infile) ~= 2
    infile = [inpath,'data',sprintf('%4.4d',vec(i)),'.asc'];
  end
  if exist(infile) ~= 2
    infile = [inpath,sprintf('%4.4d',vec(i)),'.asc'];
  end
  if exist(infile) ~= 2
    infile = [inpath,sprintf('%4.4d',vec(i)),'_data.cnv'];
  end
  if exist(infile) ~= 2
    infile = [inpath,sprintf('%4.4d',vec(i)),'_Data.cnv'];
  end
  if exist(infile) ~= 2
    infile = [inpath,sprintf('%4.4d',vec(i)),'_data.asc'];
  end
  if exist(infile) ~= 2
      ssname = regexp(moor,'_','split');
      mcatfname = [ssname{1} '_' sprintf('%4.4d',vec(i)) '_' ...
          num2str(dateoffset+1,'%4.0f') '.cnv'];
      infile = [inpath,mcatfname];
  end

  stopped=0;
  % if still doesn't exist, flag to operator
  if exist(infile) ~=2
      disp(' ')
      disp(['WARNING: No data file found for serial number ' num2str(vec(i))])
      want_to_continue=input('Do you wish to continue processing the rest of the mooring? (y/n):- ','s');
      if isempty(strmatch(upper(want_to_continue),'Y'))
          disp('Stopping!')
          stopped=1;
          fprintf(fidlog,['\n Routine stopped due to missing datafile for serial number: ' num2str(vec(i))]);
          break
      end
  else
  
    outfile = [outpath,moor,'_',sprintf('%4.4d',vec(i)),out_ext];

    % --- convert from raw to rodb format ---
    % created by Loic on PE399
    microcat2rodb_4(infile,outfile,infofile,fidlog,'y',dateoffset); 

    disp(['Press any key to continue: ']) 
    %pause
  end
end

if stopped==0
    comment = input('Enter additional comment to be save in Log file: ','s'); 

    if ~isempty(comment)
        fprintf(fidlog,'\n COMMENT:\n %s',comment);
    end
end
fclose(fidlog);
