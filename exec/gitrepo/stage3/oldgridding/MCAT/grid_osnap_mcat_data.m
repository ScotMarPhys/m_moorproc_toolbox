%=========================================================================
%=========================================================================
%%%%    GRID DATA for OSNAP mooring
%=========================================================================
%=========================================================================
close all

%pathosnap = '/home/mstar/osnap'; 
% path of the mooring data define in the startup file under osnap/
% % moor = 'rteb1_01_2014';
% moor = 'rtwb1_01_2014';
moor = 'rtwb2_01_2014';
% moor = 'rtwb1_02_2015';
% moor = 'rtwb2_02_2015';
% moor = 'rteb1_02_2015';


todo_hydrogrid       = 1;

%=========================================================================
% Grids all calibrateed moored CTD Data
% paths
p_hydrogrid.basedir      = pathosnap;
p_hydrogrid.moor         = moor;   
p_hydrogrid.ini_path     = [p_hydrogrid.basedir '/data/moor/proc/hydro_grid/ini'];
p_hydrogrid.rodbpath     = [p_hydrogrid.basedir '/data/moor/'];
p_hydrogrid.rodbctdpath  = [p_hydrogrid.basedir '/cruise_data/'];
p_hydrogrid.info_path    = [p_hydrogrid.rodbpath,'proc/',moor,'/'];
p_hydrogrid.mooringpath  = [pathosnap '/data/moor/proc'];
p_hydrogrid.out_path     = [p_hydrogrid.basedir '/data/moor/proc/hydro_grid/'];
p_hydrogrid.outname      = [moor,'_grid.mat'];
p_hydrogrid.dataexternal_ctd_dir        = [p_hydrogrid.basedir '/cruise_data/'];
p_hydrogrid.datactd_ref_cruises   = {'kn221';'pe399'};

% general settings
p_hydrogrid.dum       = -9999.0000;
p_hydrogrid.c1535     = 42.914;
p_hydrogrid.t90_68    = 1.00024;      % convert its90 to its68 for cond. to sal. conversion
p_hydrogrid.mcat      = [332:337];
p_hydrogrid.int_step  = 10;           % vertical interpolation step
p_hydrogrid.preverse  = 4000;         % 4000 pressure level below whch deep temperature reversion may occur
% repair / despike settings
p_hydrogrid.gap_max      = 10;  % allow for a maximum of gap [days] in data to interpolate across
p_hydrogrid.y_tol        = [-10 10];
p_hydrogrid.stddy_tol    = 4;
p_hydrogrid.nloop        = 5;
p_hydrogrid.graphics     = 'y';
p_hydrogrid.distrangectd = 100e3; % distance of the reference ctd from the mooring

if todo_hydrogrid  ==1
% only interpolate mooring lineary on the vertical:
%     - no use of TS climatology to get salinities where the sal. record is completely useless
%     - no guess of  the course of the profiles between the vertical grid points
hydro_grid_osnap_linear_interp(p_hydrogrid)

  if exist('/home/sa02lh')==7
    unix(['/bin/cp ' p_hydrogrid.out_path moor '_grid.mat /home/sa02lh/Data/Dropbox/Work/Postdoc_OSNAP/OSNAP_mooring/Gridded_data/.'])
    unix(['chown sa02lh:sa02lh /home/sa02lh/Data/Dropbox/Work/Postdoc_OSNAP/OSNAP_mooring/Gridded_data/'  moor '_grid.mat'])
  end
end


