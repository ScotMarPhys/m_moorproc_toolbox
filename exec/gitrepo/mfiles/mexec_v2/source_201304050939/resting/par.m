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
faco = [tstream '.ACO'];
ftpl = [tstream '.TPL'];
fmat = [tstream '.mat'];
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

num_dc_start = length(time_all);


fid = fopen(fullfaco,'r');

% find number of bytes in file
fseek(fid,0,1);
nbytes = ftell(fid);
fseek(fid,0,-1);

fseek(fid,bytes_all,-1);



block1 = 1e6; % read the file in blocks of 1e7 bytes
maxread = 1e12;

okchars = ' ,-.0123456789'; % chars expected in ACO files. tsshrp seems to have some sort of hex variable, and oceanlogger has ':'
intokchars = int8(okchars);
intcr = int8(sprintf('\r')); % char number 13
intnl = int8(sprintf('\n')); % char number 10
delim = ',';
intdelim = int8(delim);

intokchars = sort([intnl intcr intokchars]);

form = '%f,%f,%f,%f';
for kloop = 1:numdatavars
    form = [form ',%f'];
end


posnow = ftell(fid);
bytes_new = 0;


while 1
   remaining_bytes_infile = nbytes-posnow;
   remaining_bytes_requested = maxread-posnow;
   numtoread = min([block1 remaining_bytes_infile remaining_bytes_requested]);
   if numtoread < 1; break; end % no bytes left
   [sblock kount] = fread(fid,numtoread,'uint8=>int8');
   posnow = ftell(fid);
   
% % % % %    % check for unexpected chars
% % % % %    intkchars = unique(sblock);
% % % % %    intbadchars = setdiff(intkchars,intokchars);
% % % % %    
% % % % %    if ~isempty(intbadchars)
% % % % %        m1 = ['unexpected characters found in input file ' fname ' : '''];
% % % % %        m2 = char(intbadchars);
% % % % %        m3 = '''';
% % % % %        fprintf(MEXEC_A.Mfider,'%s\n',' ',[m1 m2 m3], ' ');
% % % % % %        error('error parsing ACO file');
% % % % %    end
   
   index_newlines = find(sblock == intnl); %find newlines
   numnl = length(index_newlines);
   
   if numnl < 1;break; end % no more data delimited by a newline in this fragment
   
   index_newlines = [0; index_newlines];
   data_block = nan+ones(4+numdatavars,numnl);
   data_vars = nan+ones(numdatavars,numnl);
   time_block = nan+ones(1,numnl);
   
   for kloop = 1:numnl
       sline = sblock(index_newlines(kloop)+1:index_newlines(kloop+1)-1); sline = sline(:)';
       
       % % %        sc = char(sline);
       % % %        data1 = sscanf(sc,form);
       % % %        data_block(:,kloop) = data1;

       kcr = find(sline==intcr); sline(kcr) = []; % strip out carriage return chars
       sline(sline==intcr) = []; % strip out carriage return chars
       kdelim = find(sline==intdelim);
       
       sc = char(sline);
       yyyy = str2double(sc(1:kdelim(1)-1));
       ddd = str2double(sc(kdelim(2)+1:kdelim(3)-1));
       fff = str2double(sc(kdelim(3)+1:kdelim(4)-1));
       time = datenum(yyyy,1,1)+ddd+fff-1;
       time_block(kloop) = time;

       for kloop2 = 1:numdatavars
           data_vars(kloop2,kloop) = str2double(sc(kdelim(3+kloop2)+1:kdelim(4+kloop2)-1));
       end
       


    
   end
% % %    time_block = datenum(data_block(1,:),1,1,0,0,0) + data_block(3,:) + data_block(4,:)-1;
% % %    data_vars = data_block(5:end,:);
   bytes_this_block = index_newlines(end);
   bytes_not_used = kount - bytes_this_block;
   
   time_all = [time_all time_block];
   data_all = [data_all data_vars];
   fseek(fid,-bytes_not_used,0); % back step file by the number of bytes not used
   bytes_new = bytes_new + bytes_this_block;
   bytes_all = bytes_all + bytes_this_block;
   posnow = bytes_all;
   
   save(fullfmat,'time_all','data_all','bytes_all','vnames','vunits');

   % progress report
   mean_linelength = bytes_all/length(time_all);
   bytes_remain = min([nbytes-posnow maxread-posnow]);
   approxnumremain = floor(bytes_remain/mean_linelength);
   dt = toc-toco;
   tremain = floor(approxnumremain*dt/(length(time_all)-num_dc_start));
   m1 = [sprintf('%9d',length(time_all)-num_dc_start) ' new data cycles parsed; '];
   m2 = ['approx ' sprintf('%9d',approxnumremain) ' remain.'];
   m3 = [' ' sprintf('%6d',tremain) ' seconds required.'];
   fprintf(MEXEC_A.Mfidterm,'%s\n',[m1 m2 m3]);


end

m = [sprintf('%d',posnow) ' bytes read'];
fprintf(MEXEC_A.Mfidterm,'%s\n',m)
m = [sprintf('%d',bytes_all) ' bytes parsed'];
fprintf(MEXEC_A.Mfidterm,'%s\n',m)

fclose(fid);
toc
% save(fullfmat,'time_all','data_all','bytes_all','vnames','vunits');
return

% note time for parsing one block and time for matlab i/o
