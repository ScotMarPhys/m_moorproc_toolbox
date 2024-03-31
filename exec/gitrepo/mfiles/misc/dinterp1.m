function y1 = dinterp1(x,y,x1)
[xu,iu,iv] = unique(x);
yu = y(iu);
ix = ~isnan(xu) & ~isnan(yu);
if sum(ix) > 2
    y1 = interp1(xu(ix),yu(ix),x1);
else
    y1 = NaN*x1;
end

