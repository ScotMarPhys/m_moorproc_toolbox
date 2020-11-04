function[ddmm] = dd2dm(dddd)

% Convert format from decimal degrees (dd.dd) to the degrees minutes format
% used by the Slocum Glider (ddmm.mm)
%
% eg.  -115.25 --> -11515


s = sign(dddd);      % Save original sign

dddd = abs(dddd);    % Absolute

d = fix(dddd);       % Degrees

dd = dddd - d;       % Decimal degrees

m = dd.*60;            % Decimal Degrees -> Minutes

ddmm = s.*((d*100)+m);    % Combine degrees and decimal degrees and restore sign
