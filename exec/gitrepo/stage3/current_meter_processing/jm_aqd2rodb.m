function jm_aqd2rodb(filterSpec)
% AQD2RODB Convert Nortek Aquadopp raw data to RODB format.
% Loads RAW file
% added pressure column --JM 2/2005
% allowed to load data that has no text in filename, eg. 123.dat --JM 3/2005
% changed column designations for temperature, pressure and w vel to corresspond
% to latest version of Aquadopp v1.27--JM 3/2005
% Input: CONF file (configuration file/header)
% Load: RAW file (stripped off data)
% Output: COR file (corrected for magnetic declination and sound speed) --JM 12/2005
% CORRECT FOR SOUND VELOCITY OF SEAWATER
% Aquadopp calculates for correct Sound Velocity already (get the 'sound speed used'
% column--see header)

d = dir(filterSpec);
dirname = filterSpec(1:max(findstr(filterSpec,filesep)));

for k = 1:length(d)
    confile = fullfile(dirname,d(k).name);
    header = readconf(confile);
    var = rodbhead(str2mat(header{:}));
    i = find(var == rodbhead('Filename'));
    outfile = sscanf(header{i},'%*s = %s');

    i = find(var == rodbhead('Mag_Variation'));
    magdec = sscanf(header{i},'%*s = %f');

    i = max(findstr(d(k).name,'.'));
    rawfile = sprintf('%s.raw',d(k).name(1:i-1));
    rawdat = sprintf('%s',d(k).name(1:i-1));
    rawfile = fullfile(dirname,rawfile);

    % get the data by loading it
    eval(['load ', rawfile]);
    eval(['a=X' rawdat ';']); %use for filenames that dont have text in name, eg, 123.dat, usually for Sonteks
    %eval(['a=' rawdat ';']); %use for filenames that have text eg, M123.dat
    [m,n] = size(a);
    yy=a(:,3);mm=a(:,1);dd=a(:,2);
    hh=a(:,4)+a(:,5)/60+a(:,6)/3600;

    jdate=julian(yy,mm,dd,hh);


    % CORRECT FOR SOUND VELOCITY OF SEAWATER
    % Aquadopp calculates for correct Sound Velocity already (see 'sound speed
    % used' column, see header file) using fixed S 35 and measured temperature &
    % pressure


    % CORRECT FOR MAGNETIC DECLINATION GIVEN BY MAGDEC
        uu = a(:,9)*100;  % convert m/s to cm/s
        vv = a(:,10)*100;
        t = a(:,23);
        p = a(:,21);
        %use uvrot by Visbeck
        [u,v]=uvrot_visbeck(uu,vv,magdec);
    
    % be sure to have corresponding columns for variables (u,v,t,p) when using other versions of
    % Aquadopp, I used v1.27 --JM 3/2005
    tab = [yy,mm,dd,hh,t,p,u,v];

    fmt = '%4d %2d %2d %8.4f %7.3f %7.3f %7.2f %7.2f\n';
    fid = fopen(outfile,'w');
    fprintf(fid,'%s',header{:});
    fprintf(fid,'Columns              = yy:mm:dd:hh:t:p:u:v\n');
    fprintf(fid,fmt,tab');
    fclose(fid);
end

%-------------------------------------------------------------------------------

function header = readconf(filename)
%READCONF
%  HEADER = READCONF(FILENAME)

fid = fopen(filename,'r');
i = 1;
while 1
    header{i} = fgets(fid);
    if header{i} == -1
        break
    end
    i = i + 1;
end
fclose(fid);
header = header(1:end-1);



