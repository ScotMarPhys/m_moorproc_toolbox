function [vdata vunits] = msload(instream,dn1,dn2,varlist,mode)
% function [vdata vunits] = msload(instream,dn1,dn2,mode)
%
% USE, eg
%   mtload('winch',[2009 4 4],[2009 5 4 12 0 0])
%   mtload('winch',now-0.1,now)
%   mtload('winch',now-0.1)
%   mtload winch now-0.1
%   mtload winch '[2009 4 9 12 0 0]' now
%   mtload winch '2009 4 9 12 0 0' now 'time cableout rate'
%
%   mtload(instream,dn1,dn2,varlist,mode) or
%   mtload instream dn1 dn2 varlist mode
%
% load techsas data into Matlab
% source stream is instream; can be techsas name or mexec short name.
% dn1 and dn2 are matlab datenums or datevecs that define the required 
% data period
% note silent mode; use 'q' to suppress output to the screen;
% Data are loaded from all relevant techsas files with matching stream
% name, and appended.
% varlist is a single character string and can use either variable names or
% numbers, eg '/' '1~4' 'time lat long' 'time 2~4' 'time 2 4 6'
% Use mtvars to find the variable names in a stream.
%
% 2009-09-22 fixed at noc to work with either standard (comma delimited) or 
% sed-revised (space delimited) ACO files

m_common
tstream = msresolve_stream(instream);

[mt1 mt2] = msgetdfinfo(tstream,'f'); % get time limits in case they are required for default

if ~exist('mode','var'); mode = ' '; end
if ~exist('dn1','var'); dn1 = mt1; end
if isempty(dn1); dn1 = mt1; end
if ischar(dn1); cmd =['dn1 = [' dn1 '];']; eval(cmd); end % if the arg has come in as a string, convert from char to number
if ~exist('dn2','var'); dn2 = now; end
if isempty(dn2); dn2 = now; end
if ischar(dn2); cmd =['dn2 = [' dn2 '];']; eval(cmd); end


% convert datevecs to nums; if the arguments are datenums nothing is
% changed
dn1 = datenum(dn1);
dn2 = datenum(dn2);

% get file names
fnames = mschoosefiles(tstream,dn1,dn2);
nf = length(fnames);

% get var list
[vars units] = msgetvars(tstream);
vars{end+1} = 'time'; % time is always a variable in scs
units{end+1} = 'matlab';
nv = length(vars);

% sort out the var list
if ~exist('varlist','var'); varlist = '/'; end
th.fldnam = vars;
th.noflds = nv; % create a structure equivalent to the mstar headers to parse for var names
% keyboard
varnums = m_getvlist(varlist,th);
% time always seems to be last in the techsas list; put it first if it is
% in the load list.
loadvarnames = vars(varnums);
ktime = strmatch('time',loadvarnames);
if ~isempty(ktime)
    timevarnum = varnums(ktime);
    varnums(ktime) = [];
    varnums = [timevarnum varnums];
end

% identify parts of files to load

dc1 = nan+ones(nf,1); dc2 = dc1; totdc = 0;
m = 'Counting data cycles';
if ~strncmp(mode,'q',1)
    fprintf(MEXEC_A.Mfidterm,'%s\n',m);
end

for kf = 1:nf
    fn = fnames{kf};
    fullfn = [MEXEC_G.Mscs_sed '/' fn];
    [dc1(kf) dc2(kf)] = msgetdcrange(fn,dn1,dn2);
    totdc = totdc + dc2(kf)-dc1(kf)+1;
end
       
m = [sprintf('%d',totdc) ' data cycles and ' sprintf('%d',nv) ' vars found in ' sprintf('%d',nf) ' files'];
if ~strncmp(mode,'q',1)
    fprintf(MEXEC_A.Mfidterm,'%s\n',m);
end

% now load data

for kv = varnums
    % make empty space so file doesn't grow in loop
    vuse = nan+ones(1,totdc);
    kount = 0;
    m = ['loading variable ' vars{kv}];
    if ~strncmp(mode,'q',1)
        fprintf(MEXEC_A.Mfidterm,'%s\n',m);
    end

    for kf = 1:nf
        fn = fnames{kf};
        fullfn = [MEXEC_G.Mscs_sed '/' fn];
        %     fullfn_mat = [fullfn(1:end-4) '.mat']; % replace .ACO with .mat
        %     bak for jr195: allow different read and write dirs for scs
        fullfn_mat = [MEXEC_G.Mscs_mat '/' fn(1:end-4) '.mat']; % replace .ACO with .mat
        nk = dc2(kf)-dc1(kf)+1; % load this many data cycles on this operation
%         vin = nc_varget(fullfn,vars{kv},dc1(kf)-1,nk);
        vin_cell = load(fullfn_mat);
        if strmatch(vars{kv},'time','exact')
            vin = vin_cell.time_all(dc1(kf):dc1(kf)+nk-1);
        else
            matnum = strmatch(vars{kv},vin_cell.vnames,'exact');
            if ~isempty(matnum)
                vin = vin_cell.data_all(matnum,dc1(kf):dc1(kf)+nk-1);
            else
                vin = nan+vin_cell.time_all(dc1(kf):dc1(kf)+nk-1);
            end
        end

        vuse(kount+1:kount+nk) = vin;
        kount = kount+nk;
    end
    loadvarname = vars{kv};
    loadvarname(strfind(loadvarname,'-')) = '_'; % remove minus in var name
    loadvarname(strfind(loadvarname,'/')) = '_'; % remove slash in var name
%     cmd = ['vdata.' vars{kv} ' =  vuse(:)'';']; eval(cmd);
%     cmd = ['vunits.' vars{kv} ' =  units{kv};']; eval(cmd);
    cmd = ['vdata.' loadvarname ' =  vuse(:)'';']; eval(cmd);
    cmd = ['vunits.' loadvarname ' =  units{kv};']; eval(cmd);
end

