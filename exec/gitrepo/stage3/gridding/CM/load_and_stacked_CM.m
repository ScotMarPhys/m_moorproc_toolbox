function [U,V,W,P,z] = load_and_stacked_CM(boundarydir,hydrodir,moor,JG,cm_check_plot)
    fileID = fopen([boundarydir, moor, '.dat']);
    delimiter = {'\t',' '};
    startRow = 6;
    % % Format string for each line of text:
    %   column1: text (%s)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    %   column5: double (%f)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%s%f%f%f%f%*s%*s%[^\n\r]';
    % % Read columns of data according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    file1_data = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    % % Close the text file.
    fclose(fileID);
    mooring = file1_data{1};
    sn      = file1_data{2};
    cm      = file1_data{3};
    z       = file1_data{4};
    lon     = file1_data{5};

    U = zeros(length(sn), length(JG));
    V = zeros(length(sn), length(JG));
    W = zeros(length(sn), length(JG));
    P = zeros(length(sn), length(JG));

    j=1;
    
    close all
    for i = 1: length(sn)

        infile = [hydrodir, mooring{i,:}, '_velocity_grid.mat'];
        load(infile,'dnumi','ufi','vfi','wfi','pfi');
        
        if strcmp(moor,'rteb_CM_osnap_02_2015')
            idx = find(dnumi>=datenum('22-Jun-2015 00:00:00'));
        elseif strcmp(moor,'rteb_CM_osnap_04_2017')
            idx = find(dnumi>=datenum('15-May-2017 00:00:00'));
        elseif strcmp(moor,'CM_rtwb1_osnap_02_2015')
            idx = find(dnumi>=datenum('23-Jun-2015 12:00:00'));
        elseif strcmp(moor,'CM_rtwb1_osnap_03_2016')
            idx = find(dnumi>=datenum('05-Jul-2016 12:00:00'));
        elseif strcmp(moor,'CM_rtwb1_osnap_04_2017')
            idx = find(dnumi>=datenum('12-May-2017 12:00:00'));
        elseif strcmp(moor,'CM_rtwb1_osnap_06_2020')
            idx = find(dnumi>=datenum('16-Oct-2020 00:00:00'));
        elseif strcmp(moor,'CM_rtwb2_osnap_02_2015')
            idx = find(dnumi>=datenum('23-Jun-2015 12:00:00'));
        elseif strcmp(moor,'CM_rtwb2_osnap_04_2017')
            idx = find(dnumi>=datenum('13-May-2017 12:00:00'));
        elseif strcmp(moor,'CM_rtwb2_osnap_05_2018')
            idx = find(dnumi>=datenum('11-Jul-2018 12:00:00'));
        else
            idx = 1:length(dnumi);  
        end
        
        jdnew = dnumi(idx);
        ufi =  ufi(:,idx);
        vfi =  vfi(:,idx);
        wfi =  wfi(:,idx);
        pfi =  pfi(:,idx);
        
        if strcmp(moor,'rteb_CM_osnap_03_2016')
            if cm(i)==1 || cm(i)==2
                idx = find(jdnew>=datenum('28-Mar-2017 12:00:00'));
                ufi(cm(i),idx)=NaN;
                vfi(cm(i),idx)=NaN;
                wfi(cm(i),idx)=NaN;
                pfi(cm(i),idx)=NaN;
            end
        end
        
        sampling_rate = round(1./median(diff(jdnew))); %nominal sampling rate [per day]
        uuu    = interp1(jdnew, ufi(cm(i),:), JG);
        vvv    = interp1(jdnew, vfi(cm(i),:), JG);
        www    = interp1(jdnew, wfi(cm(i),:), JG);    
        ppp    = interp1(jdnew, pfi(cm(i),:), JG);
        U(j,:) = uuu;
        V(j,:) = vvv;
        W(j,:) = www;    
        P(j,:) = ppp;

        j = j + 1;

        if cm_check_plot
            col = {'r','b','m','c','g','y','r--','b--','m--','c--','g--',...
                'y--','r:','b:','m:','c:','g:','y:','r.','b.','m.','c.',...
                'g.','y.','r'};
            figure(10011)
            hold on; box on;
            plot(JG, U(i, :), col{:,i},'LineWidth',2);
            title([moor,' - U'])

            figure(10012)
            hold on; box on;
            plot(JG, V(i, :), col{:,i},'LineWidth',2);
            title([moor,' - V'])

            figure(10013)
            hold on; box on;
            plot(JG, W(i, :), col{:,i},'LineWidth',2);
            title([moor,' - W'])

            figure(10014)
            hold on; box on;
            plot(JG, P(i, :), col{:,i},'LineWidth',2);
            title([moor,' - PRESSURE'])

        end

    end
end