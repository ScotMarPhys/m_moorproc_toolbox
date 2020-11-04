function ms_update_aco_to_mat(instream)

% convert SCS ACO file to .mat. If the file is part converted, update the conversion. 
%
% function ms_update_aco_to_mat(instream)
% 
% The mat file saves the number of bytes from the input file already
% converted.
%
% If in any doubt about the file, the best thing is to stop the
%   sed conversion from ACO in scs_raw to scs_sed
%   then delete the files for this stream in scs_sed and scs_mat
%   and then rerun the sed script and then rerun this aco_to_mat conversion
%
% In order to read data faster, all non-numeric characters are first removed by
%   sed in a conversion from scs_raw to scs_sed, eg with unix scripts
%   sedexec_stopall (to stop all sed stream conversions)
%   sedexec_startall (to restart and overwrite ACO files in scs_sed
% 
% CALLED BY: 
%   update_allmat
%   msgaps
%   mslast
%   mslistit
%   msposinfo
%   scs_to_mstar2
%
% INPUT:
%   instream: name of SCS stream to be converted or updated, in
%   cruise/data/scs_sed
%
% OUTPUT:
%   mat file for this stream, in cruise/data/scs_mat
%
% EXAMPLES:
%   ms_update_aco_to_mat('seatex-gll')
%
% UPDATED:
%   help comments added BAK 3 Jun 2014 on jr302


m_common
tstream = msresolve_stream(instream);


% tstream = 'anemometer'
% tstream = 'oceanlogger'
% tstream = 'tsshrp0'
% tstream = 'fake'
% tstream = 'gyro'
% tstream = 'ashtech'
% tstream = 'ea600'
% tstream = 'seatex-gga'
% tstream = 'seatex-gll'
% tstream = 'seatex-hdt'
% tstream = 'seatex-vtg'
% tstream = 'seatex-zda'
% tstream = 'winch'

% MEXEC_G.Mscs_root = '.';

tic; toco=toc; % start clock


% set up file names
faco = [tstream '.ACO'];
ftpl = [tstream '.TPL'];
fmat = [tstream '.mat'];
fullfaco = [MEXEC_G.Mscs_sed '/' faco];
fullftpl = [MEXEC_G.Mscs_sed '/' ftpl];
% fullfmat = [MEXEC_G.Mscs_root '/' fmat];
fullfmat = [MEXEC_G.Mscs_mat '/' fmat]; % bak for jr195: allow different read and write dirs for scs


if exist(fullfmat,'file') == 2
    load(fullfmat,'time_all','data_all','bytes_all','vnames','vunits');
    oldfile = 1;
    numdatavars = length(vnames);
else
    % load var list
    varcells = mtextdload(fullftpl);
    numdatavars = length(varcells);

    vnames = cell(numdatavars,1); % empty cells
    vunits = vnames;

    for kloop = 1:numdatavars % parse the names and units
        vcell = varcells{kloop};
        vnames{kloop} = vcell{2};
        vunits{kloop} = vcell{3};
    end
    % create empty arrays
    oldfile = 0;
    bytes_all = 0;
    data_all = [];
    time_all = [];
end

bytes_start = bytes_all;

num_dc_start = length(time_all);


fid = fopen(fullfaco,'r');

% find number of bytes in file
fseek(fid,0,1);
nbytes = ftell(fid);
fseek(fid,0,-1);

tline = fgets(fid);
com = strfind(tline,',');
fseek(fid,0,-1);

if isempty(com)
    % space delimited
    m = ['Parsing space delimited file ' fullfaco];
    fprintf(MEXEC_A.Mfidterm,'%s\n',m);
    f1 = '%f ';
else
    % comma delimited
    m = ['Parsing comma delimited file ' fullfaco];
    fprintf(MEXEC_A.Mfidterm,'%s\n',m);
    f1 = '%f,';
end


gash = fscanf(fid,f1,10*(4+numdatavars)); % read 10 lines
b10 = ftell(fid);
fseek(fid,0,-1);

mean_linelength = b10/10;
approx_lines_remain = (nbytes-bytes_all)/mean_linelength;
m1 = ['Average line length appears to be ' sprintf('%7.1f',mean_linelength) ' bytes'];
m2 = [sprintf('%d',num_dc_start) ' lines already in file'];
m3 = ['Approx ' sprintf('%d',round(approx_lines_remain)) ' lines remain to be parsed'];
fprintf(MEXEC_A.Mfidterm,'%s\n',m1,m2,m3);

kblock = 500000; % 100k lines at a time
newdata = nan + ones(numdatavars,kblock);
newtimes = nan + ones(4,kblock);


form = [f1 f1 f1 f1];
for kloop = 1:numdatavars
    form = [form f1];
end
alldone = [];

fseek(fid,0,-1);
fseek(fid,bytes_all,-1);
m = 'Starting to parse new data'; fprintf(MEXEC_A.Mfidterm,'%s\n',m);

posnow = ftell(fid);

while 1 % cycle until end of data
    tocstart = toc;
    knew = 0;
    newdata = nan+newdata;
    newtimes = nan+newtimes;
    for kline = 1:kblock % read kblock cycles at a time
%        if strcmp(tstream,'oceanlogger')  %jr239 get around as two extra columns in ocealogger
%               BAK at NOC: 17 Aug 2010
%               jr239 mod removed by BAK.
%               The way round the oceanlogger problem is to edit the TPL file in the scs_sed directory
%               The situation arises because the original ACO file has an extra time variable 
%               described in the TPL as 
%               126,oceanlogger-sampletime,YYY DDD HH:MM:SS
%               By the time this turns up in the scs_sed directory it is three variables, described by a single line of TPL.
%               So the answer is, edit the TPL so that where the original TPL has
%                  126,oceanlogger-sampletime,YYY DDD HH:MM:SS
%               the new TPL has
%                  126,oceanlogger-sampletimeyyyy,YYYY
%                  126,oceanlogger-sampletimeddd,DDD
%                  126,oceanlogger-sampletimehhmmss,HH:MM:SS
%               and this program will then correctly count the number of oceanlogger 'variables'.
%               It will also assign them sensible names and units
        [new kount] = fscanf(fid,form,4+numdatavars);
        if kount == 4+numdatavars
            knew = knew+1;
            newdata(:,knew) = new(5:end);
            newtimes(:,knew) = new(1:4);
            posnow = ftell(fid);
        else
            % end of data; truncate newdata arrays
            newdata(:,knew+1:end) = [];
            newtimes(:,knew+1:end) = [];
            fseek(fid,posnow,-1);
            alldone = 1; % set all done flag
            break
        end
    end
    % construct time, and append if more than zero new times
    yyyy = newtimes(1,:);
    ddd = newtimes(3,:);
    fff = newtimes(4,:);

    if numel(yyyy) > 0
        newtime = datenum(yyyy,1,1)+ddd+fff-1;
        time_all = [time_all newtime];
        data_all = [data_all newdata];
        m = [sprintf('%d',knew) ' new lines parsed from this block'];
%         fprintf(MEXEC_A.Mfidterm,'%s\n',m);
    end
    bytes_all = posnow; % save position after last successful data read


    m = 'Starting to save appended data';
%     fprintf(MEXEC_A.Mfidterm,'%s\n',m);toc;

    save(fullfmat,'time_all','data_all','bytes_all','vnames','vunits');
    m = 'Finished saving appended data';
%     fprintf(MEXEC_A.Mfidterm,'%s\n',m);toc;

    % progress report
    bytes_remain = nbytes-posnow;
    approxnumremain = floor(bytes_remain/mean_linelength);
    dt = toc-tocstart;
    tremain = floor(approxnumremain*dt/numel(yyyy));
    m1 = [sprintf('%9d',length(time_all)-num_dc_start) ' new data cycles parsed; '];
    m2 = ['approx ' sprintf('%9d',approxnumremain) ' remain.'];
    m3 = [' ' sprintf('%6d',tremain) ' seconds required.'];
    fprintf(MEXEC_A.Mfidterm,'%s\n',[m1 m2 m3]);

    
    
    if ~isempty(alldone); break; end

end

return



% % % % %    % progress report
% % % % %    mean_linelength = bytes_all/length(time_all);
% % % % %    bytes_remain = min([nbytes-posnow maxread-posnow]);
% % % % %    
% % % % %    dt = toc-toco;
% % % % %    tremain = floor(approxnumremain*dt/(length(time_all)-num_dc_start));


% % % % % end

m = [sprintf('%d',posnow) ' bytes read'];
fprintf(MEXEC_A.Mfidterm,'%s\n',m); 
m = [sprintf('%d',posnow-bytes_start) ' new bytes parsed'];
fprintf(MEXEC_A.Mfidterm,'%s\n',m); 

fclose(fid);
toc
% save(fullfmat,'time_all','data_all','bytes_all','vnames','vunits');
return

% note time for parsing one block and time for matlab i/o
