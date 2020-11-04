% MC_CALL_2_RB1201 is a script that performs stage1 processing
% on microcat data.  It converts microcat data from raw to rodb
% format for an entire mooring.
%
% It calls microcat2rodb_3 (to convert microcat_data), rodbload.m,
% timeaxis.m, auto_filt.m, julian.m 

% 27/10/09 DR added functionality for .cnv files
% 21/03/10 ZBS made all paths rely on a basedir variable,
%          added a descriptive header, and modified for oceanus 459
% 08/11/12 DR changed write mode of log file so appends instead of
%          overwriting
% 12/11/12 DR added check on whether infile exists and clearly highlights
%          it for the operator with the option of stopping the routine.

close all
clear all 

% -----------------------------------------------------------------
% --- This is the information that needs to be modified for -------
% --- different users, directory trees, and moorings --------------
% -----------------------------------------------------------------

%basedir  = '/Volumes/RB1201/rapid/data/'; % on RB1201
%basedir = '/noc/users/pstar/rpdmoc/rapid/data/'; % on D382
basedir = '/noc/mpoc/rpdmoc/rapid/data/'; % at NOCS after D382
cruise = 'mocha_ab1209';
operator = 'dr400';

% moor = 'wbh2_5_201116';
% moor = 'wb2_9_201114';
% moor = 'ebh5_7_201131';
% moor = 'mar2_8_201136';
%moor = 'wb6_6_201201';
%moor = 'wb1_9_201208';
%moor='wb4_9_201202';
%moor = 'wb2_10_201205';
%moor = 'wbh2_6_201204';
moor = 'mochab_6_403';
%moor = 'mochae_6_405';
%------------------------------------------


% --- set paths for data input and output ---

addpath(genpath([basedir 'exec/' cruise '/']));
%addpath(genpath([basedir 'moor']));

inpath   = [basedir 'moor/raw/' cruise '/microcat/'];
outpath  = [basedir 'moor/proc/' moor '/microcat/'];
infofile = [basedir 'moor/proc/' moor '/' moor 'info.dat'];

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
    infile = [inpath,sprintf('%4.4d',vec(i)),'_data.asc'];
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
    microcat2rodb_3(infile,outfile,infofile,fidlog,'y',0); % DR manually added 0 toffset value on D324.

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
