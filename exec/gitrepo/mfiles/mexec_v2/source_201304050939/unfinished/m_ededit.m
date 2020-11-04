% called from mplxyed

if ~exist('kfind','var')
    m = 'No data selected to edit';
    m1 = 'Select some data first with option s';
    fprintf(MEXEC_A.Mfider,'%s\n',m,m1)
    return
else
    subindex = kfind{vared}(:)';
    if length(subindex) == 0
        m = 'No data cycles selected to edit';
        m1 = 'Select some data first with option s';
        fprintf(MEXEC_A.Mfider,'%s\n',m,m1)
        return
    end
end

if kedit == 0;
    % first time in edit case
    kedit = 1; % set a flag to show that some editing has been done
    ncf = m_openio(pdfot.ncfile); % set write flag
    % create a variable to record the edits
    kedits = [];
end

index = x1-1+kfind{vared}(:)';
vnam = h.fldnam{ynumlist(vared)};
kedits = [kedits index];
kedits = unique(kedits);

vdata = nc_varget(pdfot.ncfile.name,vnam,[r1-1 c1-1],[r2-r1+1,c2-c1+1]);
vdata(kfind{vared}) = nan;
nc_varput(pdfot.ncfile.name,vnam,vdata,[r1-1 c1-1],[r2-r1+1,c2-c1+1]);
m_uprlwr(pdfot.ncfile,vnam);




%     for ki = 1:length(index)
%         vall = [];
%         for k = 1:length(vlist)
%             vnum = vlist(k);
%             vnam = h.fldnam{vnum};
%             [row col] = m_index_to_rowcol(index(ki),h,vnum);
%             vdata = nc_varget(pdfot.ncfile.name,vnam,[row-1 col-1],[1 1]);
%             vdata = reshape(vdata,numel(vdata),1);
%             %             if ~isempty(find(k_tconvert == k)) % then this var must be converted assuming it is decimal days
%             if ktime(k) == 1 & ktconvert == 1 % then this var must be converted using its units
%                 [yy mo dd hh mm ss dayofyear] = m_time_to_ymdhms(vnam,vdata,h);
%                 vdata = [yy-2000 mo dd dayofyear hh mm round(ss)]; % ss is rounded in m_time_to_ymdhms, but need to force it to integer
%                 if vdata(1) < 0; vdata(1) = vdata(1)+100; end % this will only produce a 2-digit year if 1900 <= yy <= 2099
%             end
%
%             vall = [vall vdata];
%         end
%         vall = [index(ki) vall];
%         s = sprintf(form,vall(1,:));
%         fprintf(MEXEC_A.Mfidterm,'%s\n',s);
%     end


return
