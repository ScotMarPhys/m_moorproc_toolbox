function ms_update_aco_to_mat(instream)

m_common
tstream = msresolve_stream(instream);

% matnames = msgetstreamfilenames(tstream);
% fname = matnames{1};



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
faco = [tstream 'out.ACO'];
ftpl = [tstream '.TPL'];
fmat = [tstream 'out2.mat'];
fullfaco = [MEXEC_G.Mscs_root '/' faco];
fullftpl = [MEXEC_G.Mscs_root '/' ftpl];
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
if isempty(com)
    % space delimited
    m = 'Parsing space delimited file';
    fprintf(MEXEC_A.Mfidterm,'%s\n',m);
    f1 = '%f ';
else
    % comma delimited
    m = 'Parsing comma delimited file';
    fprintf(MEXEC_A.Mfidterm,'%s\n',m);
    f1 = '%f,';
end


xx = fscanf(fid,f1,10*(4+numdatavars)); % read 10 lines
b10 = ftell(fid);
fseek(fid,0,-1);

bperline = b10/10;
approx_lines_remain = (nbytes-bytes_all)/bperline;

m1 = ['Average line length appears to be ' sprintf('%7.1f',bperline) ' bytes'];
m2 = ['Approx ' sprintf('%10.1f',approx_lines_remain) ' lines remain to be parsed'];
fprintf(MEXEC_A.Mfidterm,'%s\n',m1,m2);

lines_allow = floor(1.02*approx_lines_remain);
newdata = nan+ones(numdatavars,lines_allow);
expand = nan+ones(numdatavars,floor(0.1*approx_lines_remain));
newtimes = ones(4,lines_allow);
expand_times = nan+ones(4,floor(0.1*approx_lines_remain));


% form = '%f,%f,%f,%f,';
% for kloop = 1:numdatavars
%     form = [form '%f,'];
% end
form = [f1 f1 f1 f1];
for kloop = 1:numdatavars
    form = [form f1];
end

fseek(fid,0,-1);
fseek(fid,bytes_all,-1);
m = 'Starting to parse new data';
fprintf(MEXEC_A.Mfidterm,'%s\n',m);

knew = 0;

% suggest this is all read and saved iteratively 100k lines at a time, with
% update reports. Use sed to switch comma to space, so that space delimited
% in oceanlogger is parseable. Will need to edit TPL for oceanlogger.
% need to fix exerythign for wrong number of vars in oceanlogger and check if it works on other streams
while 1
    posnow = ftell(fid);
    [new kount] = fscanf(fid,form,4+numdatavars);
    if kount == 4+numdatavars
        knew = knew+1;
        if knew > size(newdata,2)
            newdata = [newdata expand];
            newtimes = [newtimes expand_times];
        end
        newdata(:,knew) = new(5:end);
        newtimes(:,knew) = new(1:4);
    else
        newdata(:,knew+1:end) = [];
        newtimes(:,knew+1:end) = [];
        fseek(fid,posnow,-1);
        break
    end
end
% % keyboard
% % return
m = 'Finished parsing new data';
fprintf(MEXEC_A.Mfidterm,'%s\n',m);


yyyy = newtimes(1,:);
ddd = newtimes(3,:);
fff = newtimes(4,:);

if numel(yyyy) > 0
    newtime = datenum(yyyy,1,1)+ddd+fff-1;
    time_all = [time_all newtime];
    data_all = [data_all newdata];
end
bytes_all = posnow;


m = 'Starting to save appended data';
fprintf(MEXEC_A.Mfidterm,'%s\n',m);toc;

   save(fullfmat,'time_all','data_all','bytes_all','vnames','vunits');
m = 'Finished saving appended data';
fprintf(MEXEC_A.Mfidterm,'%s\n',m);toc;

% % % % %    % progress report
% % % % %    mean_linelength = bytes_all/length(time_all);
% % % % %    bytes_remain = min([nbytes-posnow maxread-posnow]);
% % % % %    approxnumremain = floor(bytes_remain/mean_linelength);
% % % % %    dt = toc-toco;
% % % % %    tremain = floor(approxnumremain*dt/(length(time_all)-num_dc_start));
% % % % %    m1 = [sprintf('%9d',length(time_all)-num_dc_start) ' new data cycles parsed; '];
% % % % %    m2 = ['approx ' sprintf('%9d',approxnumremain) ' remain.'];
% % % % %    m3 = [' ' sprintf('%6d',tremain) ' seconds required.'];
% % % % %    fprintf(MEXEC_A.Mfidterm,'%s\n',[m1 m2 m3]);


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
