function mpaste(varargin)

% paste vars from a second file onto the first, optionally using a control variable

m_common
m_margslocal
m_varargs

MEXEC_A.Mprog = 'mpaste';
m_proghd;

fprintf(MEXEC_A.Mfidterm,'%s\n','Enter name of output disc file')
fn_ot = m_getfilename;
ncfile_ot.name = fn_ot;
ncfile_ot = m_openio(ncfile_ot);

fprintf(MEXEC_A.Mfidterm,'%s\n','Enter name of input disc file')
fn_in = m_getfilename;
ncfile_in.name = fn_in;
ncfile_in = m_openin(ncfile_in);

hin = m_read_header(ncfile_in);
hot = m_read_header(ncfile_ot);

% bug fix by bak on di346 17 feb 2010
% previously, mpaste did not put name/version of the 'output' file as one of the input files 
% in the history file.
hist = hot;
hist.filename = ncfile_ot.name;
MEXEC_A.Mhistory_in{1} = hist;

% h = m_read_header(ncfile_in);
% m_print_header(h);

hist = hin;
hist.filename = ncfile_in.name;
MEXEC_A.Mhistory_in{2} = hist;


% --------------------
% Now do something with the data
% first write the same header; this will also create the file
hot.openflag = 'W'; %ensure output file remains open for write, even though input file is 'R';
m_write_header(ncfile_ot,hot);

m1 = ['Do you want to use a control variable ?'];
m2 = ['reply ''y'' or ''n'' (default)'];
fprintf(MEXEC_A.Mfidterm,'%s\n',m1,m2);

kcontrol = 0;
okreply = 0;
while okreply == 0
    reply = m_getinput(' ','s');
    if strcmp(' ',reply) == 1; kcontrol = 0; break; end
    if strcmp('/',reply) == 1; kcontrol = 0; break; end
    if strcmp('y',reply) == 1; kcontrol = 1; break; end
    if strcmp('n',reply) == 1; kcontrol = 0; break; end
    fprintf(MEXEC_A.Mfider,'\n%s\n','You must reply one of ''y'' ''/'' return or ''n'' : ');
end

if kcontrol > 0
    ok = 0;
    while ok == 0;
        hin = m_read_header(ncfile_in);
        m_print_header(hin);
        m = sprintf('%s\n','Type variable name or number for control variable on input file for paste : ');
        var = m_getinput(m,'s');
        if strcmp(' ',var) == 1;
            vlistcin = [];
        else
            vlistcin = m_getvlist(var,hin);
        end

        if length(vlistcin) ~= 1
            m = 'You must choose precisely one control variable. try again';
            fprintf(MEXEC_A.Mfider,'%s\n',m)
            continue
        end
        vcontrolin = vlistcin;
        cdatain = nc_varget(ncfile_in.name,hin.fldnam{vcontrolin});
        cdatain = reshape(cdatain,1,numel(cdatain));
        ok = 1;
    end
    ok = 0;
    while ok == 0;
        hot = m_read_header(ncfile_ot);
        m_print_header(hot);
        m = sprintf('%s\n','Type variable name or number for control variable on output file for paste : ');
        var = m_getinput(m,'s');
        if strcmp(' ',var) == 1;
            vlistcot = [];
        else
            vlistcot = m_getvlist(var,hot);
        end

        if length(vlistcot) ~= 1
            m = 'You must choose precisely one control variable. try again';
            fprintf(MEXEC_A.Mfider,'%s\n',m)
            continue
        end
        vcontrolot = vlistcot;
        cdataot = nc_varget(ncfile_ot.name,hot.fldnam{vcontrolot});
        cdataot = reshape(cdataot,1,numel(cdataot));
        ok = 1;
    end
end

% control variable ok, now get variables for paste
listok = 0;
while listok == 0
    m_print_header(hin);
    m1 = 'Type variable names or numbers for variables from input';
    m2 = 'file for paste (return for none, ''/'' for all): ';
    m = sprintf('%s\n',m1,m2);
    var = m_getinput(m,'s');
    if strcmp(' ',var) == 1;
        vlistin = [];
    else
        vlistin = m_getvlist(var,hin);
        m = ['list is ' sprintf('%d ',vlistin) ];
        disp(m);
    end

    m_print_header(hot);
    m1 = 'Type variable names or numbers for variables from output';
    m2 = 'file for paste (return for none, ''/'' for all): ';
    m = sprintf('%s\n',m1,m2);
    var = m_getinput(m,'s');
    if strcmp(' ',var) == 1;
        vlistot = [];
    else
        vlistot = m_getvlist(var,hot);
        m = ['list is ' sprintf('%d ',vlistot) ];
        disp(m);
    end
    
    if  length(vlistin) == length(vlistot)
        % lists appear to match
        % unfinished should check dimensionsof in and out. If control var
        % is used some dimensions will need to match those if control var
        % on input file
        listok = 1;
        continue
    else
        m = 'Your lists of input and output variables do not match; try again';
        fprintf(MEXEC_A.Mfider,'%s\n',m)
    end
end

numpaste = length(vlistin);

for k = 1:numpaste
    datain = nc_varget(ncfile_in.name,hin.fldnam{vlistin(k)});
    if kcontrol == 0
        nc_varput(ncfile_ot.name,hot.fldnam{vlistot(k)},datain);
        nump = numel(datain);
        m = [sprintf('%10d',nump) ' datacycles pasted for input variable ' hin.fldnam{vlistin(k)}];
        fprintf(MEXEC_A.Mfidterm,'%s\n',m);
    else
        dataot = nc_varget(ncfile_ot.name,hot.fldnam{vlistot(k)}); % get the data before pasting
        ndata = numel(datain); % assume length of control data matches
        % unfinished, should do match once before entering variables loop
        % k. can output number of datacycles to paste (ie number matching on control var) at this stage
        nump = 0;
        for k2 = 1:ndata
            cin = cdatain(k2);
            kmat = find(cdataot == cin);
            kot = min(kmat);
            if ~isempty(kot); dataot(kot) = datain(k2); nump = nump+1; end
        end
        nc_varput(ncfile_ot.name,hot.fldnam{vlistot(k)},dataot);
        m = [sprintf('%10d',nump) ' datacycles pasted for input variable ' hin.fldnam{vlistin(k)}];
        fprintf(MEXEC_A.Mfidterm,'%s\n',m);
    end
%     varoutname = ['temp2_' sprintf('%d',k)]
    varoutname = hot.fldnam{vlistot(k)};
    %unfinished, need code to sort out output name if there is a clash, as
    %per mcalc etc
%     error(' fix code here')
%     nc_varrename(ncfile_ot.name,hot.fldnam{vlistot(k)},varoutname)
%     nc_attput(ncfile_ot.name,varoutname,'units',hin.fldunt{vlistin(k)})
% unfinished: need option of taking variable name from in or out file
% at present, var name not changed.
    m_uprlwr(ncfile_ot,varoutname);
end

%unfinished need to fix names and dimensions taken from in or out file

% % % % % %copy selected vars from the infile
% % % % % m = sprintf('%s\n','Type variable names or numbers to copy (return for none, ''/'' for all): ');
% % % % % var = m_getinput(m,'s');
% % % % % if strcmp(' ',var) == 1;
% % % % %     vlist = [];
% % % % % else
% % % % %     vlist = m_getvlist(var,h);
% % % % %     m = ['list is ' sprintf('%d ',vlist) ];
% % % % %     disp(m);
% % % % % end
% % % % % 
% % % % % for k = vlist
% % % % %     vname = h.fldnam{k};
% % % % %     numdc = h.dimrows(k)*h.dimcols(k);
% % % % %     m = ['Copying ' sprintf('%8d',numdc) ' datacycles for variable '  vname ];
% % % % %     fprintf(MEXEC_A.Mfidterm,'%s\n',m);
% % % % %     m_copy_variable(ncfile_in,vname,ncfile_ot,vname);
% % % % % end
% % % % % 
% % % % % ok = 0;
% % % % % while ok == 0;
% % % % %     m = sprintf('%s\n','Type variable name or number for control variable for merge : ');
% % % % %     var = m_getinput(m,'s');
% % % % %     if strcmp(' ',var) == 1;
% % % % %         vlistc = [];
% % % % %     else
% % % % %         vlistc = m_getvlist(var,h);
% % % % %     end
% % % % % 
% % % % %     if length(vlistc) ~= 1
% % % % %         m = 'You must choose precisely one control variable. try again';
% % % % %         fprintf(MEXEC_A.Mfider,'%s\n',m)
% % % % %         continue
% % % % %     end
% % % % %     vcontrol = vlistc;
% % % % %     ok = 1;
% % % % % end
% % % % % 
% % % % % % get control data
% % % % % vname_c_1 = h.fldnam{vcontrol};
% % % % % data_c_1 = nc_varget(ncfile_in.name,vname_c_1);
% % % % % 
% % % % % 
% % % % % m = 'Now get details of next input file ';
% % % % % fprintf(MEXEC_A.Mfidterm,'%s\n',m)
% % % % % fprintf(MEXEC_A.Mfidterm,'%s\n','Enter name of input disc file')
% % % % % fn_in2 = m_getfilename;
% % % % % ncfile_in2.name = fn_in2;
% % % % % 
% % % % % ncfile_in2 = m_openin(ncfile_in2);
% % % % % 
% % % % % h2 = m_read_header(ncfile_in2);
% % % % % m_print_header(h2);
% % % % % 
% % % % % hist = h2;
% % % % % hist.filename = ncfile_in2.name;
% % % % % MEXEC_A.Mhistory_in{2} = hist;
% % % % % 
% % % % % 
% % % % % ok = 0;
% % % % % while ok == 0;
% % % % %     m = sprintf('%s\n','Type variable name or number for control variable for merge : ');
% % % % %     var = m_getinput(m,'s');
% % % % %     if strcmp(' ',var) == 1;
% % % % %         vlistc2 = [];
% % % % %     else
% % % % %         vlistc2 = m_getvlist(var,h2);
% % % % %     end
% % % % % 
% % % % %     if length(vlistc2) ~= 1
% % % % %         m = 'You must choose precisely one control variable. try again';
% % % % %         fprintf(MEXEC_A.Mfider,'%s\n',m)
% % % % %         continue
% % % % %     end
% % % % %     vcontrol2 = vlistc2;
% % % % %     ok = 1;
% % % % % end
% % % % % 
% % % % % data_c_2 = nc_varget(ncfile_in2.name,h2.fldnam{vcontrol2});
% % % % % if m_isvartime(h2.fldnam{vcontrol2}); data_c_2 = m_adjtime(h2.fldnam{vcontrol2},data_c_2,h2,h); end % adjust time to data time origin of first input file
% % % % % 
% % % % % numdims_c = m_numdims(data_c_2);
% % % % % 
% % % % % contmerge = 0;
% % % % % 
% % % % % if numdims_c == 2
% % % % %     % 2-d control var
% % % % %     m = 'Your control var has 2 dimensions, do you want to use a row or column ';
% % % % %     m2 = 'as the independent variable for the interpolation ?';
% % % % %     fprintf(MEXEC_A.Mfidterm,'%s\n',m,m2);
% % % % %     ok = 0;
% % % % %     while ok == 0;
% % % % %         m3 = sprintf('%s','type r for row, c for column :  ');
% % % % %         reply = m_getinput(m3,'s');
% % % % %         if strcmp(reply,'r'); rc = 2; break; end % ind var is a row; colindex varies
% % % % %         if strcmp(reply,'c'); rc = 1; break; end % ind var is a col; rowindex varies
% % % % %         fprintf(MEXEC_A.Mfider,'\n%s\n','You must reply r or c : ');
% % % % %     end
% % % % % 
% % % % %     str = {'column' 'row'};
% % % % %     strtext = str{rc};
% % % % % 
% % % % %     % find row/column number
% % % % %     m1 = 'Interpolation only works when the independent variable is 1-D ?';
% % % % %     m2 = ['Which ' strtext ' do you want to use ?'];
% % % % %     fprintf(MEXEC_A.Mfidterm,'%s\n',m1,m2);
% % % % %     maxd = size(data_c_2,1);
% % % % %     ok = 0;
% % % % %     while ok == 0;
% % % % %         m3 = sprintf('%s',['type number in range 1 (default) to ' sprintf('%d',maxd) '  ']);
% % % % %         reply = m_getinput(m3,'s');
% % % % %         if strcmp(reply,' '); contindex = 1; break; end
% % % % %         cmd = ['contindex = [' reply '];']; %convert char response to number
% % % % %         eval(cmd);
% % % % %         if length(contindex) ~= 1; continue; end
% % % % %         ok = 1;
% % % % %     end
% % % % % 
% % % % %     m1 = ['Do you want to merge on the other ' strtext 's of the control variable ?'];
% % % % %     fprintf(MEXEC_A.Mfidterm,'%s\n',m1);
% % % % %     maxd = size(data_c_2,1);
% % % % %     ok = 0;
% % % % %     while ok == 0;
% % % % %         m3 = sprintf('%s',['reply no (n, default) or yes (y)  ']);
% % % % %         reply = m_getinput(m3,'s');
% % % % %         if strcmp(reply,' '); contmerge = 0; break; end
% % % % %         if strcmp(reply,'n'); contmerge = 0; break; end
% % % % %         if strcmp(reply,'y'); contmerge = 1; break; end
% % % % %     end
% % % % % 
% % % % % 
% % % % %     if rc == 1
% % % % %         x = data_c_2(:,contindex);
% % % % %     else
% % % % %         x = data_c_2(contindex,:);
% % % % %     end
% % % % % else
% % % % %     % 1-D var
% % % % %     x = data_c_2;
% % % % % end
% % % % % 
% % % % % xdm = min(diff(x));
% % % % % if xdm <= 0
% % % % %     m = 'Control variable was not monotonic; it will be sorted';
% % % % %     fprintf(MEXEC_A.Mfider,'%s\n',m)
% % % % %     [xsort ksort] = sort(x);
% % % % % else
% % % % %     xsort = x; ksort = 1:length(x);
% % % % % end
% % % % % xbad = isnan(xsort);
% % % % % if sum(xbad) > 0
% % % % %     m = ['Control variable contained '  sprintf('%d',sum(xbad)) ' nan values; these are not valid in interp1 and have been removed'];
% % % % %     fprintf(MEXEC_A.Mfider,'%s\n',m)
% % % % %     xsort(xbad) = [];
% % % % % end
% % % % % 
% % % % % % now we've got the two control variables: data_c_1 and xsort.
% % % % % % xsort has been adjusted for data time origin if necessary
% % % % % % if they're both time variables but one is days and the other is seconds,
% % % % % % scale xsort into same units as data_c_1
% % % % % 
% % % % % %vname_c_1 = vname_c_1
% % % % % unit_c_1 = h.fldunt(vlistc);
% % % % % vname_c_2 = h2.fldnam{vcontrol2};
% % % % % unit_c_2 = h2.fldunt{vcontrol2};
% % % % % 
% % % % % if (m_isvartime(vname_c_1) == 1) & (m_isvartime(vname_c_2) == 1)
% % % % %     % both recognised as time variables
% % % % %     if m_isunitsecs(unit_c_2) == 1
% % % % %         if m_isunitdays(unit_c_1) == 1
% % % % %             xsort = xsort/86400;
% % % % %         end
% % % % %     end
% % % % %     if m_isunitdays(unit_c_2) == 1
% % % % %         if m_isunitsecs(unit_c_1) == 1
% % % % %             xsort = xsort*86400;
% % % % %         end
% % % % %     end
% % % % % end
% % % % % 
% % % % % 
% % % % % % check whether data_c_1 is in range of xsort
% % % % % 
% % % % % xmin = min(xsort);
% % % % % xmax = max(xsort);
% % % % % ximin = min(min(data_c_1));
% % % % % ximax = max(max(data_c_1));
% % % % % if ximin < xmin
% % % % %     m = ['warning file 1 control variable is not contained within file 2 at low end'];
% % % % %     fprintf(MEXEC_A.Mfider,'%s\n',m)
% % % % % end
% % % % % if ximax > xmax
% % % % %     m = ['warning file 1 control variable is not contained within file 2 at high end'];
% % % % %     fprintf(MEXEC_A.Mfider,'%s\n',m)
% % % % % end
% % % % % 
% % % % % 
% % % % % % control variable ok, now get variables for merge
% % % % % 
% % % % % m = sprintf('%s\n','Type variable names or numbers for variables for merge (return for none, ''/'' for all with matching dimensions): ');
% % % % % var = m_getinput(m,'s');
% % % % % if strcmp(' ',var) == 1;
% % % % %     vlist2 = [];
% % % % % else
% % % % %     vlist2 = m_getvlist(var,h2);
% % % % %     m = ['list is ' sprintf('%d ',vlist2) ];
% % % % %     disp(m);
% % % % % end
% % % % % 
% % % % % ok = 0;
% % % % % m1 = sprintf('%s',['If NaNs are found in the merging variable, do you want them filled first (f, default)']);
% % % % % m2 = sprintf('%s',['or do you want them kept in place (k) so that NaNs may appear in the output ?']);
% % % % % fprintf(MEXEC_A.Mfidterm,'%s\n',m1,m2);
% % % % % while ok == 0;
% % % % %     m3 = sprintf('%s','reply ''f'' or return for fill or ''k'' for keep : ');
% % % % %     reply = m_getinput(m3,'s');
% % % % %     if strcmp(reply,' '); absfill = 1; break; end
% % % % %     if strcmp(reply,'f'); absfill = 1; break; end
% % % % %     if strcmp(reply,'k'); absfill = 0; break; end
% % % % % end
% % % % % 
% % % % % % find vars with 'matching' dimensions
% % % % % crows = h2.dimrows(vcontrol2);
% % % % % ccols = h2.dimcols(vcontrol2);
% % % % % ccycles = crows*ccols;
% % % % % nrows = h2.dimrows;
% % % % % ncols = h2.dimcols;
% % % % % ncycles = nrows.*ncols;
% % % % % 
% % % % % if ccols == 1; rc = 1; end; % only one col so use it as independent var;
% % % % % if crows == 1; rc = 2; end; % only one row so use it as independent var;
% % % % % 
% % % % % if(rc == 1) % we're working down columns
% % % % %     kmat = find(nrows == crows);
% % % % % end
% % % % % if(rc == 2) % we're working along rows
% % % % %     kmat = find(ncols == ccols);
% % % % % end
% % % % % 
% % % % % if contmerge == 1; vlist2 = [vcontrol2 vlist2]; end% since the user asked for the merging variable, make sure it is in the list
% % % % % 
% % % % % kmat = intersect(kmat,vlist2); % find vars that are in user's list and have suitable dims
% % % % % 
% % % % % if contmerge == 0; kmat = setdiff(kmat,vcontrol2); end % remove control var from action list
% % % % % %leave control var in if contmerge == 1
% % % % % 
% % % % % if rc == 2;
% % % % %     % need to transpose so that interp1 works on gridded data
% % % % %     xsort = xsort';
% % % % % end
% % % % % 
% % % % % m = ['list of matching vars is ' sprintf('%d ',kmat) ];
% % % % % disp(m);
% % % % % 
% % % % % 
% % % % % for k = 1:length(kmat)
% % % % %     vname = h2.fldnam{kmat(k)};
% % % % %     kmat2 = strmatch(vname,h.fldnam(vlist),'exact');
% % % % %     if ~isempty(kmat2) % attempting to copy a variable that has already been taken from first file
% % % % %         m1 = ['attempting to merge a variable                    ' vname ];
% % % % %         m2 = ['that has already been copied from first input file'];
% % % % %         m3 = ['you need to rename the variable for output'];
% % % % %         fprintf(MEXEC_A.Mfider,'%s\n',m1,m2,m3)
% % % % % 
% % % % %         ok = 0;
% % % % %         while ok == 0;
% % % % %             m3 = sprintf('%s',['type new variable name for output :              ']);
% % % % %             newname = m_getinput(m3,'s');
% % % % %             if strcmp(newname,' ') | strcmp(newname,'/');
% % % % %                 m = 'try again';
% % % % %                 fprintf(MEXEC_A.Mfider,'%s\n',m)
% % % % %                 continue
% % % % %             end
% % % % %             newname = m_remove_outside_spaces(newname);
% % % % %             newname = m_check_nc_varname(newname);
% % % % %             ok = 1;
% % % % %         end
% % % % %     else
% % % % %         newname = vname;
% % % % %     end
% % % % % 
% % % % %     m = ['Merging ' h2.fldnam{kmat(k)}];
% % % % %     fprintf(MEXEC_A.Mfidterm,'%s\n',m);
% % % % %     z = nc_varget(ncfile_in2.name,h2.fldnam{kmat(k)});
% % % % %     if m_isvartime(h2.fldnam{kmat(k)}); z = m_adjtime(h2.fldnam{kmat(k)},z,h2,h); end % adjust time to data time origin of first input file
% % % % %     % rearrange data according to sort of control variable
% % % % %     if rc == 2;
% % % % %         % need to transpose so that interp1 works on gridded data
% % % % %         z = z';
% % % % %     end
% % % % %     zsort = z(ksort,:);
% % % % %     zsort(xbad,:) = [];
% % % % %     
% % % % %     if absfill == 1;
% % % % %         % fill any nans in z by interpolation
% % % % %         for kfill = 1:size(zsort,2);
% % % % %             ok = ~isnan(zsort(:,kfill));
% % % % %             if sum(ok) < 2 ; continue; end % not enough good data to fill with interp1
% % % % %             if sum(ok) == size(zsort,1); continue; end % all data are good; skip interp1
% % % % %             zsort(:,kfill) = interp1(xsort(ok),zsort(ok,kfill),xsort);
% % % % %         end
% % % % %     end
% % % % %     
% % % % %     zi = interp1(xsort,zsort,data_c_1); % zi has same dimensions as  data_c_1
% % % % %     % yi = interp1(x,y,xi)
% % % % %     % if y is 2-D then yi is always a column, so we need to fix its shape
% % % % %     % to match xi
% % % % %     % if y is 1-D, then yi is the same shape as xi
% % % % %     if m_numdims(zsort) == 1 % zsort is 1-D
% % % % %         ktranspose = 0; % dimensions of zi will match data_c_1;
% % % % %     else % zsort is 2-D and zi will have merge index varying down a column 
% % % % %         if size(data_c_1,1) == 1
% % % % %             % data_c_1 was a row vector. Force z to be a row vector
% % % % %             ktranspose = 1;
% % % % %         else
% % % % %             ktranspose = 0;
% % % % %         end
% % % % %     end
% % % % %     
% % % % %     if ktranspose == 1
% % % % %         zi = zi';
% % % % %     end
% % % % % 
% % % % % 
% % % % %     % write the output
% % % % %     clear v
% % % % %     v.name = newname;
% % % % %     v.data = zi;
% % % % %     m_write_variable(ncfile_ot,v,'nodata'); %write the variable information into the header but not the data
% % % % %     % next copy the attributes
% % % % %     vinfo = nc_getvarinfo(ncfile_in2.name,vname);
% % % % %     va = vinfo.Attribute;
% % % % %     for k2 = 1:length(va)
% % % % %         vanam = va(k2).Name;
% % % % %         vaval = va(k2).Value;
% % % % %         nc_attput(ncfile_ot.name,newname,vanam,vaval);
% % % % %     end
% % % % % 
% % % % %     % now write the data, using the attributes already saved in the output file
% % % % %     % this provides the opportunity to change attributes if required, eg fillvalue
% % % % %     nc_varput(ncfile_ot.name,newname,v.data);
% % % % %     m_uprlwr(ncfile_ot,newname);
% % % % % 
% % % % % end







% finish up

m_finis(ncfile_ot);

h = m_read_header(ncfile_ot);
m_print_header(h);

hist = h;
hist.filename = ncfile_ot.name;
MEXEC_A.Mhistory_ot{1} = hist;
m_write_history;


return






