% script to cd to the directory stored in global variable MEXEC_G.MEXEC_CWD
% set MEXEC_G.MEXEC_CWD with function 'mcsetd'
%
% USE:
%  mcd;
%
% INPUT:
%   none; MEXEC_G.MEXEC_CWD must be set with mcsetd.
%
% OUTPUT:
%   none; matlab current working directory is changed
%
% UPDATED:
%   Initial version BAK 2008-10-17 at NOC


m_common
MEXEC_A.Mprog = 'mcd';
m_proghd;


clear swhos
swhos = whos('MEXEC_G.MEXEC_CWD');

if isempty(swhos)
    fprintf(MEXEC_A.Mfider,'%s\n','MEXEC_G.MEXEC_CWD not set. Set it using ''mcsetd'' '); return
elseif swhos.global == 0
         fprintf(MEXEC_A.Mfider,'%s\n','MEXEC_G.MEXEC_CWD does not appear to be global. Set it using ''mcsetd'' '); return
elseif exist(MEXEC_G.MEXEC_CWD,'dir') ~= 7
    fprintf(MEXEC_A.Mfider,'%s\n',['Directory selected by MEXEC_G.MEXEC_CWD does not exist: ' MEXEC_G.MEXEC_CWD]); return
else
    cmd = ['cd ' MEXEC_G.MEXEC_CWD ';']; eval(cmd)
    m = ['Directory changed to ' MEXEC_G.MEXEC_CWD];
    fprintf(MEXEC_A.Mfidterm,'%s\n',m);
end



