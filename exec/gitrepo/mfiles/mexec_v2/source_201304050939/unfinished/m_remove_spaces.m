function sot = m_remove_spaces(sin)

% remove any spaces from a string

s = sin;
k = strfind(sin,' ');
s(k) = [];
sot = s;
return