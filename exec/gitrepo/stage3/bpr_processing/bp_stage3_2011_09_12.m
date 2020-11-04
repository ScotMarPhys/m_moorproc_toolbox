% BP_STAGE3_2010-03-01 is a script to process BP records that need
% stage 3 processing as of 2010-03-01.  Stage 3 processing performs
% harmonic fitting and subtracts weekly and monthly tides.
%

% SKYE 2012/02/24; cleaned up a bit
% Z Szuts 2011/08/31
% Z Szuts 2010/08/03
% Z Szuts 2010/01/27
% Z Szuts 31/1/09


% operator specific paths etc

basedir = '~/Desktop/RB1201/rapid/data/moor/proc/';
setenv('ROBDPATH',basedir')

% moorings to process for RB1201,
%moors = {...};
moors = {'wb6_5_201117'};

outvars = ['Latitude:Longitude:WaterDepth:StartDate:StartTime:'...
           'EndDate:EndTime:z:instrument:serialnumber'];


% valid instrument types for BPR sensors (for seagauge, ixsbpr,
% wlr, seacat, bourdon, and pies)
bpr_ids = [465 470 460 332 480 316]; 


for i=1:length(moors)
  moor = moors{i};

  disp([' --- mooring ' moor ' --- '])
  disp(' ')

  infofile = [basedir moor '/' moor 'info.dat'];

  disp(['infofile: ', infofile]);
  [lat,lon,wd,sdate,stime,edate,etime,z,type,sn] = ...
      rodbload(infofile,outvars);


  isbp = zeros(size(type));
  for j = 1:length(type)
    if any( type(j) == bpr_ids );
      isbp(j) = 1;
    end
  end

  js = find(isbp);

  % only process one type of sensor at a time, but
  % bottom_pressure_grid can process multiple sensors of the same
  % type with one call
  if length(js) ~= length(unique(type(js)))
    jjs = js(1);
    types = type(js(1));

    js(1) = [];
    while length(js)
      if ~ismember(type(js(1)),types)
        jjs = [jjs js(1)];
      end
      js(1) = [];
    end      
    js = jjs;
    clear types jjs
  end    

  
  for j = js

    if type(j) == 465
      ext = '.seagauge';
    elseif type(j) == 470
      ext = '.ixsbpr';
    elseif type(j) == 460
      ext = '.wlr';
    elseif type(j) == 332
      ext = '.seacat';
    elseif type(j) == 480
      ext = '.bourdon';   
    elseif type(j) == 316
      ext = '.pies';
    else
      error('the value of type is not an expected BPR instrument')
    end  
    bp_id = type(j);
    
    % just call bottom_pressure_grid to generate stage 3 files
    %[p_grid,jd_grid,tide_fit] = bottom_pressure_grid(moor,basedir,type(j));

    bottom_pressure_grid3 % now in script form
    clear moor bp_id
  end

  disp(' ')
  disp(' ')

end

