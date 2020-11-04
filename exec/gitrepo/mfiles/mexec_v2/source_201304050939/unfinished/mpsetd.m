function mpsetd(varargin)
% function mpsetd(varargin)
%
% can be used in a matlab session or called from a script
%
% similar to mcsetd; discovers the full path to a directory used for
% datapup reading from 'rvs' files to 'pstar' files

m_common
m_margslocal
m_varargs

MEXEC_A.Mprog = 'mpsetd';
m_proghd;

% m_getinput takes input from MEXEC_A.MARGS_IN via varargin if present

m1 = 'Set a directory as the pexec current working directory';
m2 = ['This will be below the dataroot being used for datapup work: ' MEXEC_G.Mrsh_dataroot_local];
m = sprintf('%s\n',m1,m2);
dn = m_getinput(m,'s');

% sort out if it is one in the global list

mstar_list = MEXEC_G.PDIRLIST(:,1);
target_list = MEXEC_G.PDIRLIST(:,2);

k = strmatch(dn,mstar_list,'exact');
if isempty(k)
    % not one of the global list, assume it is just a simple directory name
else
    dn = target_list{k};
end


MEXEC_G.PEXEC_CWD_FULL = [MEXEC_G.Mrsh_dataroot_local '/' dn]; % set current working directory; set relative path name only
MEXEC_G.PEXEC_CWD = [dn]; % set current working directory; set relative path name only
m2 = ['MEXEC_G.PEXEC_CWD set to : ' MEXEC_G.PEXEC_CWD];
fprintf(MEXEC_A.Mfidterm,'%s\n',m2);

if exist(MEXEC_G.PEXEC_CWD_FULL,'dir') ~= 7
    fprintf(MEXEC_A.Mfider,'%s\n',['Warning: Directory selected by MEXEC_G.PEXEC_CWD does not exist: ' MEXEC_G.PEXEC_CWD_FULL]);
end

return

