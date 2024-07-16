function[d,m] = dd2dm2(dddd)

% Convert format from decimal degrees (dd.dd) to 
% degrees and minutes as two variabels.
%
% eg.  -115.25 --> -115 15.00


s = sign(dddd);      % Save original sign
dddd = abs(dddd);    % Absolute
d = fix(dddd);       % Degrees

dd = dddd - d;       % Decimal degrees

m = dd.*60;            % Decimal Degrees -> Minutes
d = s*d;
