function mhistory

% type the history entry for the most program

m_common


cmd = ['!tail -'  sprintf('%d',MEXEC_A.Mhistory_lastlines) ' ' MEXEC_A.Mhistory_filename]
eval(cmd)