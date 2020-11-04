% PDRIFT_COMPARE_OC459 is a script to compare the records of BPRs
% recovered on OC459 with previous records, to see whether the
% initial exponential transient is reproduced (and thus is oceanic)
% in independent measurements or not.

% 03 Apr 2010, ZB Szuts, oc459


% -----------------------------------------------------------------
% --- This is the information that needs to be modified for -------
% --- different users, directory trees, and moorings --------------
% -----------------------------------------------------------------

cruise   = 'oc459';
operator = 'zszuts';
mooring1  = 'wbl3_2_200806';  sn1 = 28;
mooring2  = 'wbl1_2_200705';  sn2 = [398 399];

plot_interval = [2007  1  1 00;   % start time of time axis on plot
                 2010 12 31 00];  % end time of time axis on plot

% oceanus 459
% determine whether this script is being run on Brian King's
% computer (dhcp108), or else remotely by a computer that has
% mounted the surman directory.
[out,host] = unix('hostname');
if strfind(host,'dhcp108.science.oceanus.whoi.edu');
  basedir  = '/Users/surman/rpdmoc/rapid/data/';
else
  basedir  = '/Volumes/surman/rpdmoc/rapid/data/';
end

% -----------------------------------------------------------------


% --- set paths for data input and output ---
% NB The seagauge dir in outpath must be created first

inpath1 = [basedir 'moor/proc/' mooring1 '/seagauge/'];
inpath2 = [basedir 'moor/proc/' mooring2 '/seagauge/'];
outpath = inpath1;

% initialize figure
figure
hold on, box on, grid on
set(0,'defaulttextinterpreter','none')

dum = -9999; % dummy value

jd0 = julian(-1,1,1,24);
jd1 = julian(plot_interval(1,:))-jd0;
jd2 = julian(plot_interval(2,:))-jd0; 

% set colors for all BPR records
cols1 = [0 0 0; % black
         0.5 0.5 0.5]; % gray
cols2 = [0 0 1; % blue
         0 1 0]; % green

legendtext = {};

for i = 1:length(sn1)
  
  mooringtext = [mooring1 '_' sprintf('%4.4d',sn1(i))];
  infile1 = [inpath1 mooringtext '.use'];

  [YY,MM,DD,HH,P,Pfit,T] = rodbload(infile1,'YY:MM:DD:HH:P:PFIT:T');
  jd = julian(YY,MM,DD,HH);
  
  sampling_rate = round(1./median(diff(jd))); % nominal sampling rate [per day] 
  
  Pf = repmat(nan,size(P));
  ii = find(isfinite(P) & P~=dum);
  Pf = auto_filt(P(ii),sampling_rate,1/2,'low',4);

  plot(jd-jd0,Pf-Pfit,'-','Color',cols1(i,:))
  
  disp(['Look at the stage2 logfile for ' mooringtext ':'])
  pfittext = [];
  while isempty(pfittext)
    pfitopt = input('   Is a (1) linear+exponential or a (2) linear fit used? ');
    if pfitopt==1
      pfittext = 'lin+exp';
    elseif pfitopt==2
      pfittext = 'lin';
    else
      disp('invalid input, enter again')
    end
  end
  
  legendtext{end+1} = [mooringtext ', ' pfittext];
end


for i = 1:length(sn2)
  
  mooringtext = [mooring2 '_' sprintf('%4.4d',sn2(i))];
  infile2 = [inpath2 mooringtext '.use'];

  [YY,MM,DD,HH,P,Pfit,T] = rodbload(infile2,'YY:MM:DD:HH:P:PFIT:T');
  jd = julian(YY,MM,DD,HH);
  
  Pf = repmat(nan,size(P));
  ii = find(isfinite(P) & P~=dum);
  Pf = auto_filt(P(ii),sampling_rate,1/2,'low',4);
  
  plot(jd-jd0,Pf-Pfit,'-','Color',cols2(i,:))

  disp(['Look at the stage2 logfile for ' mooringtext ':'])
  pfittext = [];
  while isempty(pfittext)
    pfitopt = input('   Is a (1) linear+exponential or a (2) linear fit used? ');
    if pfitopt==1
      pfittext = 'lin+exp';
    elseif pfitopt==2
      pfittext = 'lin';
    else
      disp('invalid input, enter again')
    end
  end

  legendtext{end+1} = [mooringtext ', ' pfittext];

end

hl = legend(legendtext);
set(hl,'Interpreter','none')

%xlim([0 julian(plot_interval(2,:))-jd0])
%datetick('x',12)

xlim([jd1 jd2])
datetick('x',12)


ylabel('Pressure (db)')
title(['Comparing BPR de-drifted data from: ' mooring1 ' and ' ...
       mooring2])

orient landscape



outfile = [outpath mooring1 '_Pdrift_compare'];
print(gcf,'-depsc',[outfile '.eps']) 
