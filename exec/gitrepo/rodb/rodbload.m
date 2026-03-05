function varargout = rodbload(ident,variables)
%RODBLOAD Load RO database files.
%  [V1,V2,...,VN] = RODBLOAD(IDENT,VARIABLES)
%
% input  :      ident           : string containing filenames
%                                 either 'clear' or 'coded' 
%               variables       : string containing up to 15 variable
%                                 identifiers
%
% output :      v1            : first variable
%               ...
%               vn            : n-th variable
%
%------------------------------------------------------------------------------
%
% TYPE 'rodbhelp' FOR MORE INFORMATIONS ON ARGUMENTS
% TYPE 'roident' FOR MORE INFORMATIONS ON IDENTIFIERS
%
% uses :        rodbhead.m, str2pos.m
%
% version 0.6.0         last change 26.02.1999
%
%
% C. Mertens and G. Krahmann, IfM Kiel, Sep 1995
%
% Update LOH, 04/12/2017, to add compatbility for windows system path (see
% line 115)

% initial coding                 C. Mertens, Sep 1995,   0.0.1-->0.2.6
% new keyword recognition        G.Krahmann, Oct 1995,   0.2.6-->0.3.0
% removed bugs                   G.Krahmann, Oct 1995,   0.3.1-->0.3.2
% removed bug                    G.Krahmann, Jan 1996,   0.3.2-->0.3.3
% added multiple search paths    C. Mertens, Jan 1996,   0.3.3-->0.3.4
% cleaned up the code            C. Mertens, Jan 1996,   0.3.4-->0.3.5
% further clean up               G.Krahmann, Jan 1996,   0.3.5-->0.3.6
% removed bug                    G.Krahmann, Apr 1996,   0.3.6-->0.3.7
% added RODBPATH                 C. Mertens, Apr 1996,   0.3.7-->0.3.8
% minor bug removed              C. Mertens, Apr 1996,   0.3.8-->0.3.9
% bug removed                    C. Mertens, Apr 1996,   0.3.9-->0.4.0
% minor bug removed              C. Mertens, Apr 1996,   0.4.0-->0.4.1
% allowed for `.dat'             C. Mertens, Apr 1996,   0.4.1-->0.4.2
% fifteen output variables       C. Mertens, Jun 1996,   0.4.2-->0.4.3
% support for multisegment files C. Mertens, Jan 1997,   0.4.3-->0.4.9
% allowed for `.raw'             C. Mertens, May 1997,   0.4.9-->0.5.1
% don't abort on emtpy file name C. Mertens, May 1997,   0.5.1-->0.5.2
% bug removed                    C. Mertens, May 1997,   0.5.2-->0.5.3
% updated to Matlab 5            C. Mertens, Feb 1999,   0.5.3-->0.6.0


global allkeys allkeys1 nall


%-------------------------------------------------------------------------------
%
% initialization
%
%-------------------------------------------------------------------------------

varargout = cell(1,nargout);

% delimiter for tokens in `ident' and `variables'
d = ':,';

% root directory of data base, multiple paths could be specified by setting
% the environment variable `RODBPATH' or by using `rodbpath'
dbpath = getenv('RODBPATH');
if exist('rodbpath') == 2
  dbpath = rodbpath;
end
if length(dbpath) == 0
  dbpath = ''; %dbpath = '/home/hamlet/ro2';  % modified by loic houpert July 2016
end

% a single file name could be given in `ident'
if length(findstr(ident(1,:),d(1))) < 2
  onefile = 1;
  ndir = 1;
else
  onefile = 0;
  ndir = size(ident,1);
end

% get keyword numbers of wanted informations
  [varnumber,varformat,varlength] = rodbhead(variables);


% initialize output arrays
for i = 1:length(varnumber)
  eval(['v',int2str(i),' = [];'])
end

% prepare maximum data-length vector
rowmax = ones(1,length(varnumber));

colnum = rodbhead('columns:columns cont.');


%-------------------------------------------------------------------------------
%
% loop over all directories
%
%-------------------------------------------------------------------------------

% counter for columns
j = 0;

for k = 1:ndir

  if onefile
    nfiles = size(ident,1);
  else
      if ispc
          % if windows distribution add an additional split of the filepath
          % due to the colon sperator for drive symbol (e.g. M:\)
             % get directory of experiment
            [driveletter,str] = strtok(ident(k,:),d);    
             % get directory of experiment
            [experiment,str] = strtok(str,d);      
      else
            % get directory of experiment
            [experiment,str] = strtok(ident(k,:),d);          
      end
    % get cruise or other subdirectory
    [cruise,str] = strtok(str,d);
    % get instrument
    [instrument,str] = strtok(str,d);
    % get file identifiers
    eval(['file = ',str(2:length(str)),';']);
    % number of files to load
    nfiles = length(file);
  end

  % prepare output arrays
  for i = 1:length(varnumber)
    if varnumber(i) >= 0
      m = max(size(varargout{i},1),1000);
      varargout{i} = [varargout{i} NaN*ones(m,nfiles)];
    end
  end


%-------------------------------------------------------------------------------
%
% loop over all files in the current experiment and instrument
%
%-------------------------------------------------------------------------------

  ifile = 0;
  while nfiles > 0

    ifile = ifile + 1;
    nfiles = nfiles - 1;

    % get file name
    if onefile
      filename = deblank(ident(ifile,:));
    else
      % set name of directory
      dirname = [experiment,'/',cruise,'/',instrument,'/'];
      % compose filename
      filename1 = [dirname,cruise,'_',sprintf('%3.3d',file(ifile)),'.', ...
			  instrument];
      
      filename2 = [dirname,cruise,'_',sprintf('%3.3d',file(ifile)),'.dat'];
      filename3 = [dirname,cruise,'_',sprintf('%4.4d',file(ifile)),'.use'];
      filename4 = [dirname,cruise,'_',sprintf('%3.3d',file(ifile)),'.raw'];
      filename5 = [dirname,cruise,'_',sprintf('%3.3d',file(ifile)),'.cal'];
      filename6 = [dirname,cruise,sprintf('%3.3d',file(ifile)),'.', ...
			  instrument];
      
      % search path for data file
      pathtmp = dbpath;
      filename = [];
%       while length(pathtmp) > 0
%         [rootdir,pathtmp] = strtok(pathtmp,d);
%         if exist([rootdir,'/',filename1]) == 2
%           filename = [rootdir,'/',filename1];
%         elseif exist([rootdir,'/',filename2]) == 2
%           filename = [rootdir,'/',filename2];
%         elseif exist([rootdir,'/',filename3]) == 2
%           filename = [rootdir,'/',filename3];
%         elseif exist([rootdir,'/',filename4]) == 2
%           filename = [rootdir,'/',filename4];
%         elseif exist([rootdir,'/',filename5]) == 2
%           filename = [rootdir,'/',filename5]
%         elseif exist([rootdir,'/',filename6]) == 2
%           filename = [rootdir,'/',filename6];
%         end
% modified by Loic Houpert (July 2016), seems to be an old definition of the path
% according to Kiel filesystem:
        if exist(filename1) == 2
          filename = filename1;
        elseif exist(filename2) == 2
          filename = filename2;
        elseif exist(filename3) == 2
          filename = filename3;
        elseif exist(filename4) == 2
          filename = filename4;
        elseif exist(filename5) == 2
          filename = filename5
        elseif exist(filename6) == 2
          filename = filename6;
        end
% keyboard
    end

    % open input file
    if isempty(filename)
      fp = -1;
      message = 'Empty file name.';
    else
      fprintf('%s\n',filename)
      [fp,message] = fopen(filename,'r');
    end
    if fp == -1
      message = sprintf('rodbload: %s\n%s\n',filename,message);
      fprintf('%s\n',message)
    else


%-------------------------------------------------------------------------------
%
% read header lines; lines with a `#' or 'remarks' at the beginning of a line
% are interpreted as comments
%
%-------------------------------------------------------------------------------

      cont = 1;
      while cont

        % initialize some counters
        count = 0;
        columns = '';
        allheaders = [];

        while count == 0
          % read header line
          headerline = fgetl(fp);
          if strcmp(headerline(1:3),'tru') | strcmp(headerline(1:3),'trD')
              headerline = ['Ins',headerline];
          end 
         % store position indicator
          fpos = ftell(fp);
          % store header lines
        
          allheaders = str2mat(allheaders,headerline);

       % read data; if this fails continue with header lines
          [data,count] = fscanf(fp,'%f');
        end

        % remove first empty line from headerline array
        allheaders = allheaders(2:size(allheaders,1),:);
  
        % reread first data line and get number of data columns
        fnext = ftell(fp);
        if feof(fp)
          cont = 0;
        end
        fseek(fp,fpos,'bof');
        line = fgetl(fp);
        fseek(fp,fnext,'bof');
        ncols = length(sscanf(line,'%f'));
        m = length(data);

        % reshape data
        data = reshape(data,ncols,m/ncols)';

        % commented out on RB0701 as causing problems with replacing valid
        % currents data. NaNs then not subsequently able to be loaded with
        % rodbload.!!!
        % replace dummy values with NaN's
       % i = find(abs(data - 1e32) < eps);
       % ij = find(abs(data - [-9.0]) < eps);
       % ik = find(abs(data - [-9.9]) < eps);
       % data(i) = NaN*data(i);
       % data(ij) = NaN*data(ij);
       % data(ik) = NaN*data(ik);


%------------------------------------------------------------------------
%
% rearrange `data' to match identifiers
%
%------------------------------------------------------------------------

        % get keyword numbers of header lines
        headnumber = rodbhead(allheaders);
        % get keyword numbers of data
        colheads = find(headnumber == colnum(1) | headnumber == colnum(2));
        datanumber = [];
        for i = 1:length(colheads)
          j1 = find(allheaders(colheads(i),:) == '=');
          col = allheaders(colheads(i),j1+1:size(allheaders,2));
          j1 = min(find(abs(col) ~= 32));
          j2 = max(find(abs(col) ~= 32));
          col = col(j1:j2);
          datanumber = [datanumber,rodbhead(col)];
        end

        % find out matching numbers and give data
        j = j + 1;
        for i = 1:length(varnumber)
          if varnumber(i) ~= -9999
            ind = find(varnumber(i) == headnumber);
            if ~isempty(ind)
              [dummy,headinfo] = strtok(allheaders(ind,:),'=');
              mi = min(find(abs(headinfo(2:length(headinfo))) ~= 32));
              if isempty(mi)
                mi = 2;
              end
              headinfo = headinfo(mi:length(headinfo));

              if strcmp(deblank(varformat(i,:)),'pos')
                varargout{i}(j) = str2pos(headinfo);
              else
                if varlength(i) == 1
                  varargout{i}(j) = sscanf(headinfo,varformat(i,:));
                else
                  dummy = sscanf(headinfo,varformat(i,:));
                  if isstr(dummy)
                    varargout{i}{j} = dummy;
                  else
                    varargout{i}(1:varlength(i),j) = dummy;
                  end
                end
              end
            else
              ind = find(varnumber(i) == datanumber);
              if ~isempty(ind)
                nlines = size(data,1);
                nv = size(varargout{i},1);
                if nv < nlines
                  varargout{i}(nv+1:nlines,1:j) = NaN*ones(nlines-nv,j);
                else
                  varargout{i}(nlines+1:nv,j) = NaN*ones(nv-nlines,1);
                end
                varargout{i}(1:nlines,j) = data(:,ind);
                rowmax(i) = max([rowmax(i),nlines]);
              end
            end
          end
        end

      end

      % close data file
      fclose(fp);

    end
  end

end

% remove spare NaNs
for i = 1:length(varnumber)
  if varnumber(i) >= 0
    varargout{i} = varargout{i}(1:rowmax(i),1:j);
  end
end
