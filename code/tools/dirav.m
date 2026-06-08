function [mang] = dirav(dang)
ang = dang*pi/180;
sang = sin(ang);
cang = cos(ang);
msang = nanmean(sang);
mcang = nanmean(cang);
mang = atan2(msang,mcang)*180/pi;
