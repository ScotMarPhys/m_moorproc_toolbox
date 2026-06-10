function datt = rodb2masked_table(datc, rodb_varstr, dummy)
% follows call to rodbload:
% datc = rodbload(infile, rodb_varstr);
% turn (single, cell) output of rodbload into a table, with NaNs in
% non-time variables

vars = strsplit(rodb_varstr,':');
msk = ~ismember(vars,{'yy','mm','dd','hh'});
datt = table(datc{1},'VariableNames',vars(1));
for vno = 1:length(vars)
    d = datc{vno};
    if msk(vno)
        d(d==dummy) = NaN;
    end
    datt.(vars{vno}) = d;
end
datt.jd = julian(datt.yy, datt.mm, datt.dd, datt.hh);
