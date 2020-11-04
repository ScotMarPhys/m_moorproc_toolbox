function fn = m_getfilename(ncfile)

% get the name of an mstar file. If there is an argument, unpack it; if
% not, prompt the keyboard

m_common

if nargin == 1 % take filename from argument
    fname = m_resolve_filename(ncfile);
    fn = fname.name;
% % % % %     fname = ncfile;
% % % % %     if isstruct(fname)
% % % % %         % assume the real file name is field 'name' of struct variable ncfile
% % % % %         fn = fname.name;
% % % % %     else
% % % % %         %its just a simple file name
% % % % %         fn = fname;
% % % % %     end
else % prompt from keyboard
    msg = 'Type name of mstar file ';
    fn = m_getinput(msg,'s','no_ot');
end


fn = m_add_nc(fn);
MEXEC_A.MARGS_OT = [MEXEC_A.MARGS_OT fn];

return