function medita(varargin)

% Edit data that lie outside a range to absent value
% varlist can be names or variable numbers

% function medita(ncfile,varlist,rangelist)

% example:
% medita(ncfile,'time u v',[0 -100 -200; 1000 100 200])
% or just type medita and be prompted for input

prog = 'medita';
m_proghd(prog);
m_global_args;


fn = m_getfilename(varargin{:}); % this inserts the optional argument if there is one
ncfile.name = fn;
ncfile = m_ismstar(ncfile);

ncfile = m_openio(ncfile);

% ncfile.name = m_add_nc(ncfile.name);
h = mlisth(ncfile);


if nargin == 3
    %     check number of vars in list matches number of range pairs
    %     varlist must be 1 x N structure array
    %     rangelist must be 2 x N double array
    vlist = m_getvlist(varlist,h);
    numedit = length(vlist);
    if size(rangelist,1) ~= 2
        error('if submitting the range list as an argument its dimension must be 2 x N');
    end

    if size(rangelist,2) ~= numedit
        errstr1 = 'if submitting the range list as an argument its dimension must be 2 x N';
        errstr2 = 'where N is the number of variables to edit';
        errstr3 = ['list of vars to edit has length ' sprintf('%2d',numedit) ' and is: ' sprintf('%d ',vlist)];
        errstr4 = ['dimensions of limits array are ' sprintf('%3d x ',size(rangelist,1)) sprintf('%5d',size(rangelist,2))];
        errstr5 = sprintf('\n%s\n%s\n%s\n%s\n',errstr1,errstr2,errstr3,errstr4);
        error(errstr5)
    end

    % now edit all vars together

    for k = 1:numedit
        vdata = nc_varget(ncfile.name,h.fldnam{vlist(k)});
        lim = rangelist(:,k);

        kbad = find(vdata < lim(1) | vdata > lim(2));
        disp(['Variable ' sprintf('%-10s',h.fldnam{vlist(k)}) ' Editing ' sprintf('%9d',length(kbad)) ' data cycles outside range ' sprintf('%12.4f %12.4f',lim)]);
        vdata(kbad) = nan;

        nc_varput(ncfile.name,h.fldnam{vlist(k)},vdata);

        m_uprlwr(ncfile,h.fldnam{vlist(k)});

    end


elseif nargin ~= 1
    error('Must call medita with precisely 1 or 3 arguments')
else
    %prompt for each variable and limits
    while 1 > 0
        varlist = input('Type variable name or number (return to finish): ','s');
        if isempty(varlist); break; end
        vlist = m_getvlist(varlist,h);
        low = input(['Type lower limit for variable ' sprintf('%d',vlist') ' (return for -inf) : ']);
        if isempty(low); low = -inf; end
        upr = input(['Type upper limit for variable ' sprintf('%d',vlist') ' (return for +inf) : ']);
        if isempty(upr); upr = inf; end

        vdata = nc_varget(ncfile.name,h.fldnam{vlist});
        lim = [low upr];

        kbad = find(vdata < lim(1) | vdata > lim(2));
        disp(['Variable ' sprintf('%-10s',h.fldnam{vlist}) ' Editing ' sprintf('%9d',length(kbad)) ' data cycles outside range ' sprintf('%12.4f %12.4f',lim)]);
        ok = input('OK to proceed (c/r for yes, anything else for no) ? ','s');
        if isempty(ok)
            vdata(kbad) = nan;
        else
            disp('                                                          Skipping edit')
        end

        nc_varput(ncfile.name,h.fldnam{vlist},vdata);

        m_uprlwr(ncfile,h.fldnam{vlist});
        h = m_read_header(ncfile);

        m_print_varsummary(h);

    end

end


m_finis(ncfile);

h = mlisth(ncfile);

return