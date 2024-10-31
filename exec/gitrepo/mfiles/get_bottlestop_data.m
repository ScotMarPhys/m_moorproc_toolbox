function [mc, ctd] = get_bottlestop_data(bottle, ctd, mc, pd, calp, ninst)
%
% extract MC and CTD data at bottle stops (stored in bottle)
% calp contains tolerances


%identify single value for each bottle stop
bst2 = [1; find(diff(bottle.p)<-calp.bottlestop_dpmin)+1];
bot_start = bottle.datnum(bst2);
nstop = length(bot_start);

if strcmp(calp.apply_offset,'y')
  bot_start = bot_start - calp.offset;
  bot_end   = bot_end - calp.offset;
  fprintf(1,' O F F S E T   H A S   B E E N   A P P L I E D  ! ! !\n')
elseif  strcmp(calp.apply_offset,'n')
  fprintf(1,' O F F S E T   H A S   N O T   B E E N   A P P L I E D  ! ! !\n')
elseif  strcmp(calp.apply_offset,'i')
  fprintf(1,' D A N G E R !!! INDIV.  OFFSETS  HAVE  BEEN   APPLIED  ! ! !\n')
  mc.datnum  = mc.datnum + ones(size(mc.datnum,1),1)*calp.impact_offset' ;
end


% ---- MicroCAT data ------------------------

calp.interval_move1=calp.interval_move/3600/24;

mcatbotstart0=nan(nstop,ninst);
mcatbotend0=nan(nstop,ninst);
mcatbotstart=nan(nstop,ninst);
mcatbotend=nan(nstop,ninst);
for stop = 1 : nstop    % bottle_stops loop
    for inst = 1 : ninst    % instrument loop

        dp=gradient(mc.P(:,inst));
        dcond = gradient(mc.C(:,inst));

        % find nearest time in MC record to bottle stop satrt time
        [~,indstop] = nearest(bot_start(stop),mc.datnum(:,inst));
        % match that time with a pressure
        presstop = mc.P(indstop,inst);
        imcatbotok00 = find(mc.P(:,inst)>presstop-3 & mc.P(:,inst)<presstop+3 ...
            & mc.datnum(:,inst)>bot_start(stop)+calp.interval_move1(1) ...
            & mc.datnum(:,inst)<bot_start(stop)+calp.interval_move1(2));% & abs(dcond)<0.02 );
        if ~isempty(imcatbotok00)
            % Add a condition to remove the first  30sec of the bottle stop.
            mdntime0 = mc.datnum(imcatbotok00(1),inst);
            imcatbotok = imcatbotok00(mc.datnum(imcatbotok00,inst)>mdntime0 + 0.5/24/60);
            % And check that the length of the bottlestop is at least == to calp.bottlestop_tmin (in sec)

            if ~isempty(imcatbotok) && ((mc.datnum(imcatbotok(end),inst) - mc.datnum(imcatbotok(1),inst))*3600*24>calp.bottlestop_tmin)

                mcatbotstart0(stop,inst)=mc.datnum(imcatbotok(1),inst);
                mcatbotend0(stop,inst)=mc.datnum(imcatbotok(end),inst);
            end
        end
    end

    diffmcatbot = mcatbotend0(stop,:)-mcatbotstart0(stop,:);
    imcsel = find(diffmcatbot==min(diffmcatbot),1);

    for inst = 1 : ninst
        if isempty(mcatbotstart0(stop,imcsel))
            mcatbotstart(stop,inst) = nan;
            mcatbotend(stop,inst) = nan;
        else
            mcatbotstart(stop,inst) = mcatbotstart0(stop,imcsel);
            mcatbotend(stop,inst)   = mcatbotend0(stop,imcsel);
        end
    end
end

figure; plot((mc.datnum(:,:)-mc.datnum(1,1))*24*60,mc.P(:,:))

for stop = 1 : nstop    % bottle_stops loop
  
  figure(10+stop); clf; hold on; ooo = 200/60/60/24;

   for inst = 1 : ninst   

       ii_move = find(mc.datnum(:,inst)<=mcatbotend(stop,inst) & ...
                      mc.datnum(:,inst)>=mcatbotstart(stop,inst));
       ii_move2 = find(mc.datnum(:,inst)<=(mcatbotend(stop,inst)+2*ooo) & ...
                      mc.datnum(:,inst)>=(mcatbotstart(stop,inst)-2*ooo));                 
       mc.tav(stop,inst) = mean(mc.T(ii_move,inst));
       mc.cav(stop,inst) = mean(mc.C(ii_move,inst));
       mc.pav(stop,inst) = mean(mc.P(ii_move,inst)); 
       if calp.cstat == 1 
         plot((mc.datnum(ii_move2,inst)-mc.datnum(1,1))*24*60,mc.C(ii_move2,inst),'k') 
         plot((mc.datnum(ii_move,inst)-mc.datnum(1,1))*24*60,mc.C(ii_move,inst),'b') 
       else
         plot((mc.datnum(ii_move2,inst)-mc.datnum(1,1))*24*60,mc.P(ii_move2,inst),'k') 
         plot((mc.datnum(ii_move,inst)-mc.datnum(1,1))*24*60,mc.P(ii_move,inst),'b')          
       end    
   end
   title([calp.cruise,'   calp.cast',num2str(calp.cast),'  depth: ',num2str(round(bottle.p(stop)))])
end

% ------  extract bottlestop values from ctd 

if rms(calp.interval_move) ~= 0
  if  calp.apply_offset   == 'y'
    ctd.time =  - offset + ctd.datnum; 
  else 
    ctd.time =  ctd.datnum;
  end
end

if calp.ctd_latestart_offset ~=0
    ctd.time = ctd.time + calp.ctd_latestart_offset/86400;
    fprintf(1,'CTD LATESTART OFFSET HAS BEEN APPLIED !!!\n')
end

ctdbotstart = nan(nstop,1);
ctdbotend = nan(nstop,1);

for stop = 1 : nstop

    ctdbotstart(stop) = min(mcatbotstart(stop,:));
    ctdbotend(stop) = max(mcatbotend(stop,:));

    ii_move = find(ctd.time<=ctdbotend(stop)...
        & ctd.time>=ctdbotstart(stop));
    ii_move2 = find(ctd.time<=ctdbotend(stop)+ooo...
        & ctd.time>= ctdbotstart(stop)-2*ooo);
    figure(10+stop)
    if calp.cstat == 1
        if strcmp('d334',calp.cruise)
            plot((ctd.time(ii_move2)-mc.datnum(1,1))*24*60,ctd.cond(ii_move2)*10,'r')
            plot((ctd.time(ii_move)-mc.datnum(1,1))*24*60,ctd.cond(ii_move)*10,'b')
        elseif isfield(pd,'ctd1hz_cunit')==1 && strcmp(pd.ctd1hz_cunit,'S/m')
            plot((ctd.time(ii_move2)-mc.datnum(1,1))*24*60,ctd.cond(ii_move2)*10,'r')
            plot((ctd.time(ii_move)-mc.datnum(1,1))*24*60,ctd.cond(ii_move)*10,'b')
        else
            plot((ctd.time(ii_move2)-mc.datnum(1,1))*24*60,ctd.cond(ii_move2),'r')
            plot((ctd.time(ii_move)-mc.datnum(1,1))*24*60,ctd.cond(ii_move),'b')
        end
    else
        plot((ctd.time(ii_move2)-mc.datnum(1,1))*24*60,ctd.temp(ii_move2),'r')
        plot((ctd.time(ii_move)-mc.datnum(1,1))*24*60,ctd.temp(ii_move),'b')
    end
    ctd.pav(stop) = mean(ctd.press(ii_move));
    ctd.tav(stop) = mean(ctd.temp(ii_move));
    if strcmp('d334',calp.cruise)
        ctd.cav(stop) = mean(ctd.cond(ii_move)*10); % pstar units are S/m and not mS/cm like normal
    elseif exist('ctd_1hz')==1 && strcmp(ctd_1hz,'S/m')
        ctd.cav(stop) = mean(ctd.cond(ii_move)*10); % cases where 1hz data in S/m and not mS/cm
    else
        ctd.cav(stop) = mean(ctd.cond(ii_move));
    end
    grid on
end
