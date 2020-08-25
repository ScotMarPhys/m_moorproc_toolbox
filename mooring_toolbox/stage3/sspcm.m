function SSP=sspcm(PTS,T,S)
% function ssp=sspcm(p,t,s)
%
% sound speed in sea water after Chen and Millero
%
% input  :	p		: pressure [dbar]
%		t		: in situ temperature [degrees Celsius]
%		s		: salinity [psu]
%
% output :	ssp		: sound speed [m/s]
%
% version 1.0.0		last change 01.09.1995
%
% for compatibility reasons it is possible to give [p,t,s] instead of
% p as only 1 argument !

%	CHECK-VALUE:1731.9954 M/S FOR P=1000 BARS,T=40 DEG C,S=40 NSSU
%	IFM KIEL	T.MUELLER
%	matlab conversion  Uwe Send, IfM/Kiel Jan. 1991
%	pts-options  23/02/93, Gerd Krahmann
%	revised header		G.Krahmann, IfM Kiel, Sep 1995
	
	if (nargin==1)
  	  P=PTS(:,1); T=PTS(:,2); S=PTS(:,3);
	else
	  P=PTS;
	end

	P=P/10.;
	SR=sqrt(abs(S));
%  S**2 TERM
	D=1.727E-3-7.9836E-6*P;
%	S**3/2 TERM
	B1=7.3637E-5+1.7945E-7*T;
	B0=-1.922E-2-4.42E-5*T;
	B=B0+B1.*P;
%  S**1 TERM
	A3=(-3.389E-13*T+6.649E-12).*T+1.100E-10;
	A2=((7.988E-12*T-1.6002E-10).*T+9.1041E-9).*T-3.9064E-7;
	A1=(((-2.0122E-10*T+1.0507E-8).*T-6.4885E-8).*T-1.2580E-5).*T+9.4742E-5;
	A0=(((-3.21E-8*T+2.006E-6).*T+7.164E-5).*T-1.262E-2).*T+1.389;
	A=((A3.*P+A2).*P+A1).*P+A0;
%  S**0 TERM
	C3=(-2.3643E-12*T+3.8504E-10).*T-9.7729E-9;
	C2=(((1.0405E-12*T-2.5335E-10).*T+2.5974E-8).*T-1.7107E-6).*T+3.1260E-5;
	C1=(((-6.1185E-10*T+1.3621E-7).*T-8.1788E-6).*T+6.8982E-4).*T+0.153563;
	C0=((((3.1464E-9*T-1.47800E-6).*T+3.3420E-4).*T-5.80852E-2).*T+5.03711).*T+1402.388;
	C=((C3.*P+C2).*P+C1).*P+C0;
%  SOUND SPEED RETURN
	SSP=C+(A+B.*SR+D.*S).*S;
