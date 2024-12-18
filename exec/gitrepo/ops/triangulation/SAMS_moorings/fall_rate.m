% fall_rate.m

load rtadcp_fall_rate.txt

yy=rtadcp_fall_rate(:,1);
mn=rtadcp_fall_rate(:,2);
dd=rtadcp_fall_rate(:,3);
hh=rtadcp_fall_rate(:,4);
mm=rtadcp_fall_rate(:,5);
ss=rtadcp_fall_rate(:,6);


range=rtadcp_fall_rate(:,7);


dnum = datenum([yy,mn,dd,hh,mm,ss]) ;
clear dt dr fr;
dt=gradient(dnum);
dt=dt.*(24*60);
dr=gradient(range);dt=gradient(dnum);
fr=dr./dt;

figure(1);clf;hold on;grid on;
plot(dnum,range,'ko-')
datetick('x')
set(gca,'YDIR','reverse')

figure(2);clf;hold on;grid on;
plot(dnum,fr,'ko-')
datetick('x')
set(gca,'YDIR','reverse')
