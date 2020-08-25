function D = dsal(XR,XT)
%DSAL	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM, Kiel

	D=((((13.5405*XR-28.1044).*XR+42.2823).*XR+50.7702).*XR ...
                  -0.1692)+(XT./(1.0+0.0162*XT)).*((((-0.0720*XR+0.2544).*XR ...
                  -0.1125).*XR-0.0132).*XR-0.0056) ;
