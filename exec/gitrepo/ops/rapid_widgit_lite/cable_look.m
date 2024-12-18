function cable_look

m_setup
t1 = datenum([2018 10 22 0 0 0]);
t2 = datenum([2017 10 22 0 0 0]); 
t3 = datenum([2017 10 22 0 0 0]); 

usb = mtload('usbpos',now-0.75,now);
win = mtload('winch' ,now-0.75,now);
t_usb = usb.time  + MEXEC_G.Mtechsas_torg-10/86400;
t_win = win.time  + MEXEC_G.Mtechsas_torg-10/86400;

dt = 300/86400;
tmx = [t1:dt:t_usb(end)];

for i = 1:length(tmx)-1
   ixa = find(t_usb > tmx(i) & t_usb < tmx(i+1));
   alt(i) = nanmedian(usb.alt(ixa));
   ixb = find(t_win > tmx(i) & t_win < tmx(i+1));
   cab(i) = nanmedian(win.cablout(ixb));
   ten(i) = nanmedian(win.tension(ixb));
end

tmy = tmx(2:end)-dt/2;
cab_rate = diff(cab)./(24*diff(tmy));
alt_rate = diff(alt)./(24*diff(tmy)); 

figure
subplot(3,2,1)
plot(tmy(2:end),cab_rate,'r')
hold on;grid on
plot(tmy(2:end),-alt_rate,'b')
ylim([-800 0])
xlim([t1 t3])

legend('Cable','Beacon','Location','North')
xlim([t1 t3])
set(gca,'XTick',[t1:1/4:t3])
datetick('x',15,'keeplimits','keepticks')
ylabel('Rate (m/hr)')
xlabel('Time (UTC)')

subplot(3,2,2)
plot(tmy,cab,'r')
hold on;grid on
xlim([t1 t3])
set(gca,'XTick',[t1:1/4:t3])
datetick('x',15,'keeplimits','keepticks')
ylabel('Cable out (m)')
xlabel('Time (UTC)')


rel_dat = [
2017 3 8 02 58 0 4937 
2017 3 8 08 28 0 3461 
2017 3 8 08 44 0 3350 
2017 3 8 09 11 0 3273 
2017 3 8 09 39 0 3228
2017 3 8 10 30 0 3173
2017 3 8 11 21 0 3100
2017 3 8 12 17 0 3006 
2017 3 8 13 47 0 2795 
2017 3 8 14 44 0 2647
2017 3 8 15 38 0 2422
2017 3 8 16 33 0 1947
2017 3 8 16 55 30 1728
2017 3 8 17 20 30 1408 
2017 3 8 17 51 00 1182
2017 3 8 18 16 00 1065
2017 3 8 19 25 00 979
];
tvec = rel_dat(:,1:6);
t_rel = datenum(tvec);
rel_dep = rel_dat(:,7);

alt2 = dinterp1(tmy,-alt,t_rel);

subplot(3,2,3)
plot(tmy,ten,'r')
hold on;grid on
xlim([t1 t3])
set(gca,'XTick',[t1:1/4:t3])
datetick('x',15,'keeplimits','keepticks')
ylabel('Cable tension (T)')
xlabel('Time (UTC)')


subplot(3,2,5)
plot(t_rel,rel_dep-alt2,'b-+')
grid on
set(gca,'YDir','reverse')
xlim([t1 t3])
set(gca,'XTick',[t1:1/4:t3])
datetick('x',15,'keeplimits','keepticks')
ylabel('Dpth CTD rel. beac (m)')
xlabel('Time (UTC)')
 
subplot(3,2,[4 6])
plot(cab,ten,'r')
hold on;grid on
ylabel('Cable tension (T)')
xlabel('Cable out (m)')
xlim([0 5000])
ylim([0 7])  