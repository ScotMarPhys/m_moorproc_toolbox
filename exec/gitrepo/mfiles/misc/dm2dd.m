function[dddd] = dm2dd(ddmm)
% Convert degrees inutes to decimal degrees
s = sign(ddmm);         % Save original sign
ddmm = abs(ddmm);       % Absolute
d = fix(ddmm/100);      % Degrees
m = ddmm-(d*100);       % Minutes
dd = m/60;              % Minutes -> Decimal Degrees
dddd = s.*(d+dd);       % Combine degrees and decimal degrees and restore sign
