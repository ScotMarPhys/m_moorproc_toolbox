%
% SCU JC064 : Copy of lines from microcat2_rodb_3 for converting time in a
% cnv file to a date and time. 
    toffset = 0;
    
    cnv_secs = input('Enter time in seconds past 2000: ', 's');
if isempty(cnv_secs)
    cnv_secs = '346786201';
end
cnv_secs = str2num(cnv_secs);
    jd=cnv_secs/(60*60*24)+julian(2000,1,1,0) - toffset
    gtime=gregorian(jd)
    HH=hms2h(gtime(:,4),gtime(:,5),gtime(:,6))
    Start_Date = gtime(1,[1 2 3])
%     End_Date = gtime(dtl,[1 2 3])