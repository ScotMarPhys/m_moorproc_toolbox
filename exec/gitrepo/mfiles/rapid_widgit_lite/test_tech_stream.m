function [tts] = test_tech_stream(stream,vars)
	% Check data stream
m_common

switch MEXEC_G.Mshipdatasystem
case 'techsas'
	strnm = mtresolve_stream(stream);
case 'rvdas'
	strnm = mrresolve_table(stream);
case 'scs'
	strnm = msresolve_stream(stream);
end

% If can't find then this will be true
if strcmp(strnm,stream)
	fprintf(1,' **** unknown stream %s ***** \n',stream)
else
	switch MEXEC_G.Mshipdatasystem
	case 'techsas'
		stvars = mtgetvars(stream);
	case 'rvdas'
		stvars = mrgetvars(stream);
	case 'scs'
		stvars = msgetvars(stream);
	end
	fprintf(1,'Steram vars are:       %s \n',char(stvars)')
	fprintf(1,'Looking for variables: %s\n',vars)

	switch MEXEC_G.Mshipdatasystem
	case 'techsas'
			[d, ~] = mtload(stream,now-0.002,now,vars);
	case 'rvdas'
			[d, ~, ~] = mrload(stream,now-0.002,now,vars);
	case 'scs'
			[d, ~] = msload(stream,now-0.002,now,vars);
	end
	if isempty(d)
		fprintf(1,' **** No data for stream %s **** \n ',stream) 
		tts = 1;
	else
		fprintf(1,' Stream %s ok \n ',stream) 
		tts = 0;
	end
end
