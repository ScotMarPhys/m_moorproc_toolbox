function S = sal(XR,XT)
%SAL	Required by SAL78.

%	10-Mar-93, C. Mertens, IfM Kiel

	S = ((((2.7081*XR-7.0261).*XR+14.0941).*XR+25.3851).*XR ...
     	    -0.1692).*XR+0.0080+(XT./(1.0+0.0162*XT)).*(((((-0.0144*XR+ ...
     	     0.0636).*XR-0.0375).*XR-0.0066).*XR-0.0056).*XR+0.0005) ;

