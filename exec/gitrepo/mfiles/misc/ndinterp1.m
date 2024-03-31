function idata = ndinterp1(zdata,yvar,tmx)
  if ~isempty(zdata)
      zdat = zdata.(yvar);
      if sum(~isnan(zdat)) > 2
         idata =  dinterp1(zdata.time,zdat,tmx);  
      else
         idata = NaN * ones(size(tmx));
      end
  else
    idata = NaN * ones(size(tmx));
  end
  