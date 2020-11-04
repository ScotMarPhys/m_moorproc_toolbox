function ms_generate_varname_translation(instream)

% generate a template file for translating variable names from an SCS stream
%
% function ms_generate_varname_translation(instream)
% SCS variable names in the TPL file include the name of the SCS variable
% This leading string is removed, eg
% seatex_gll_lon
% becomes
% lon
%
% INPUT:
%   instream: name of SCS input stream
%
% OUTPUT:
%   no output arguments
%   a rename file is created in cruise/data/templates directory,
%   eg scs_jr302_renamelist_seatex_gll.csv
%
% EXAMPLES:
%   ms_generate_varname_translation('seatex-gll')
%
% UPDATED:
%   help comments added BAK 3 Jun 2014 on jr302

m_common
tstream = msresolve_stream(instream);
tunder = tstream;
tunder(strfind(tunder,'-')) = '_';
tunder(strfind(tunder,'.')) = '_';

% set up file names
ftpl = [tstream '.TPL'];
fullftpl = [MEXEC_G.Mscs_sed '/' ftpl];

if exist(fullftpl,'file') ~= 2; return; end % skip if no TPL file

mcsetd('M_TEMPLATES','q'); root_template = MEXEC_G.MEXEC_CWD;
fntemplate = [root_template '/' MEXEC_G.Mshipdatasystem '_' MEXEC_G.MSCRIPT_CRUISE_STRING '_renamelist_' tunder '.csv'];

% load var list
varcells = mtextdload(fullftpl);
numdatavars = length(varcells);

vnames = cell(numdatavars,1); % empty cells
vunits = vnames;

m = ['Generating var name translation table for stream ' tstream];
fprintf(MEXEC_A.Mfidterm,'%s\n',m);

fid = fopen(fntemplate,'w');

for kloop = 1:numdatavars % parse the names and units
    vcell = varcells{kloop};
    vnames{kloop} = vcell{2};
    vunits{kloop} = vcell{3};

    lenstream = length(tstream);
    oldname = vnames{kloop};
    oldname(strfind(oldname,'-')) = '_';
    oldname(strfind(oldname,'/')) = '_'; % bak on jr281 there is a slash in the emlog f/a name
    newname = oldname(lenstream+2:end);
    if strcmp(newname,'time')
        % there is an instrument variable name called time that will clash
        % with the SCS time variable, so use a different name
        newname = 'inst_time';
    end
    newunits = vunits{kloop};
    if strcmp(newunits,' '); newunits = 'number'; end
    fprintf(fid,'%s,%s,%s\n',oldname,newname,newunits);
end

fclose(fid);
