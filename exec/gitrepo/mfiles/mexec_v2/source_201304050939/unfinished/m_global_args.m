% m_global_args
% 
% mexec script
% declares MEXEC_A to be global
%
% sets fields
% % output fidterm = 1 is screen; fider = 2 is standard error printed in red;
% MEXEC_A.Mfidterm = 1;
% MEXEC_A.Mfider = 2;
% % identify which variable names will be recognised as time
% MEXEC_A.Mtimnames = {'time'}; % any variable whose name begins 'time' eg 'time' 'timenew'
% identify time units that will be recognised as days and assumed to be relative to
% mstar_data_origin
% MEXEC_A.Mtimunits_days = {'day'}; % any unit whose name begins 'day' eg 'day' 'day_of_year'
% identify time units recognised as seconds and assumed to be relative to
% mstar_data_origin
% MEXEC_A.Mtimunits_seconds = {'sec'}; % any unit whose name begins 'sec' eg 'sec' 'seconds' 'sec_of_year'

global MEXEC_A

% global MEXEC_A.MARGS_IN MEXEC_A.MARGS_OT MEXEC_A.MARGS_IN_OLD MEXEC_A.MARGS_IN_LOCAL MEXEC_A.MARGS_IN_LOCAL_OLD
% global MEXEC_A.Mhistory_in MEXEC_A.Mhistory_ot MEXEC_A.Mprog MEXEC_A.Mhistory_filename MEXEC_A.Mhistory_lastlines
% global MEXEC_A.Mhistory_skip
% global MEXEC_A.Mfidterm MEXEC_A.Mfider
% global MEXEC_A.MSAVE_VLIST
% global MEXEC_A.Mtimnames MEXEC_A.Mtimunits_days MEXEC_A.Mtimunits_seconds

% output fidterm = 1 is screen; fider = 2 is standard error printed in red;
MEXEC_A.Mfidterm = 1;
MEXEC_A.Mfider = 2;

% variable names recognised as time
% MEXEC_A.Mtimnames = {'time' 'tim' 'other_string'}
MEXEC_A.Mtimnames = {'time'}; % any variable whose name begins 'time' eg 'time' 'timenew'
% time units recognised as days and assumed to be relative to
% mstar_data_origin
MEXEC_A.Mtimunits_days = {'day'}; % any unit whose name begins 'day' eg 'day' 'day_of_year'
% time units recognised as seconds and assumed to be relative to
% mstar_data_origin
MEXEC_A.Mtimunits_seconds = {'sec'}; % any unit whose name begins 'sec' eg 'sec' 'seconds' 'sec_of_year'