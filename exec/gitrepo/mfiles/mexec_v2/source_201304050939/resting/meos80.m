function mskeleton_mcalc(varargin)


% script to convert u,v to speed,direction
% using call to mcalc

m_global_args;

MEXEC_A.Mprog = 'mskeleton_mcalc';
m_proghd;
MEXEC_A.MARGS_IN = [varargin MEXEC_A.MARGS_IN];

prog = MEXEC_A.Mprog; % save for later
% output fidterm = 1 is screen; fider = 2 is standard error printed in red;
fidterm = 1;
fider = 2;


fprintf(1,'%s\n','Enter name of input disc file')
fn_in = m_getfilename;
ncfile_in.name = fn_in;
ncfile_in = m_openin(ncfile_in);
h = m_read_header(ncfile_in);
m_print_header(h);
hist = h;
hist.filename = ncfile_in.name;
history_in{1} = hist;


fprintf(1,'%s\n','Enter name of output disc file')
fn_ot = m_getfilename;
ncfile_ot.name = fn_ot;

m = sprintf('%s\n','Type variable names or numbers to copy (return for none, ''/'' for all): ');
varcopy = m_getinput(m,'s');


%--------------------------------


m1 = ' Which option ? ';
m2 = ' 1 ptmp';
m3 = ' 2 Convert Speed & Direction to E & N';
m = sprintf('%s\n',' ',m1,m2,m3);
var = m_getinput(m,'s');
if strcmp(var,' ') == 1; var = '1'; end
if strcmp(var,'1') == 1


m = 'Type variable names or numbers of east & north speed: ';
m1 = sprintf('%s\n',m);
var2 = m_getinput(m1,'s');
enlist = m_getvlist(var2,h);

m = 'speed variable name : ';
name1 = m_getinput(m,'s');
m = 'speed units (return to use same as ''east'') : ';
units1 = m_getinput(m,'s');
if strncmp(' ',units1,1) == 1; units1 = h.fldunt{enlist(1)}; end

eq1 = 'y = sqrt(x1.*x1 + x2.*x2)';

%--------------------------------
MEXEC_A.MARGS_IN = {
    fn_in
    fn_ot
    varcopy
    var2
    eq1
    name1
    units1
    '0'
    };

margs = MEXEC_A.MARGS_OT; % keep a record of input arguments for this prog
MEXEC_A.Mhistory_skip = 1; % don't write a history from the call to mcalc
% mcalc(fn_in,fn_ot,varcopy,var2,eq1,name1,units1,var2,eq2,name2,units2,'0') % can also call mcalc with arguments
mcalc
MEXEC_A.Mhistory_skip = 0;
%--------------------------------

MEXEC_A.Mprog = prog;
MEXEC_A.Mhistory_in = history_in; % retrieve value for this prog
MEXEC_A.MARGS_OT = margs;
h = m_read_header(ncfile_ot);
% m_print_header(h);

hist = h;
hist.filename = ncfile_ot.name;
MEXEC_A.Mhistory_ot{1} = hist;
m_write_history;
