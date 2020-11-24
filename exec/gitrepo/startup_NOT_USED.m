%% startup file for cruise data processing

%  CRUISE DY120 2017

%%% adds paths to subroutines and data files
% efw: updated from dy053 to dy078 on 04/05/2017
% gdm: updated from rb1201 to jc103 on 23/4/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp ('Welcome back to Matlab!')
disp ('-----------------------------------------------------------------------------------------------')
disp('    ')
disp ('DY120, 2020')
disp ('-----------------------------------------------------------------------------------------------')



disp (['this is the dy120 startup file to open the OSNAP paths to '...
    'osnap/data/moor/ and osnap/exec/dy120, '])


basedir = [ pathosnap filesep];  %'/home/mstar/osnap/'; 

%m_setup
%this_cruise = MEXEC_G.MSCRIPT_CRUISE_STRING;
this_cruise = 'dy120';
baselogdir = [basedir];
addpath(genpath([basedir 'exec/' this_cruise '/']))
%addpath(genpath([basedir 'data/moor/proc/']));
%addpath(genpath([basedir 'data/moor/proc_calib/' this_cruise '/']));
%addpath(genpath([basedir 'data/moor/raw/' this_cruise '/']));


% addpath([basedir 'amoc/m/psor/']);
% addpath([basedir 'amoc/m/tsor/']);
% addpath([basedir 'amoc/grdat/']);
% addpath([basedir 'amoc/grout/']);
% addpath([basedir 'amoc/grmoc/']);


if 0
    % setup a diary to log the matlab sessions
    diaryname = [baselogdir filesep 'matlab_diaries' filesep 'diary-' ...
        datestr(now,['yyyymmddHHMMSS'])];
    diary(diaryname);

    disp('this session is being logged to:')
    disp(diaryname)
    disp(sprintf('\n\n\n'))
end
