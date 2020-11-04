function [afit,yfit] = exp_quad_fit_seq3(x,y,aest,int_lin,int_exp)
% sequential exponential quadratic fit:
%   yfit = a1*exp( (x-x(1))/a2 ) + a3*(x-mean(x)) + a4
% of some function y = y(x) using a sequential method designed for
% fitting microcat pressure sensors.  The coefficients mean:
%   a1 = the exponential amplitude (positive for positive pressure
%        perturbations)
%   a2 = the exponential decay period (negative for decay)
%   a3 = the linear slope with the exponential removed
%   a4 = the value of the quadratic fit at mean(x)
% N.B.: the coefficients are defined differently than in exp_lin_fit2.m.
%
% The multi-step fitting procedure is:
%  1) remove a quadratic fit over int_lin (latter portion of data)
%  2) fit an exponential over int_exp (initial portion of data)
%  3) remove the exponential and fit a second quadratic over entire
%     record
%  4) fit a second exponential over int_exp
%  5) fit a third quadratic to the entire record
%  6) total fit is the second exponential plus the third quadratic fit
%
% returns coefficients of model a = [a(1) a(2) a(3) a(4)];
%
% function [a,yfit] = exp_quad_fit_seq3(x,y[,aest,int_lin,int_exp])
%
%   x       -- independent variable (units of days is assumed)
%   y       -- data points to fit
%   aest    -- initial guess for coefficients of fit,
%              a = [a1 a2 a3 a4], default = [0 0 0 0]
%   int_lin -- the portion over which to fit the initial
%              quadratic, default is [x(1)+100 x(end)]
%   int_exp -- the portion over which to fit the exponential decay,
%              default default is [x(1) x(1)+100]
%
%   afit    -- a = [a1 a2 a3 a4]: coefficients of fit 
%   yfit    -- yfit = a1*exp( (x-x(1))/a2 ) + a3*(x-mean(x)) + a4
%
%  for more help, see build-in function 'nlinfit.m'
%
%
% based partly on exp_lin_fit2.m
%
% Z Szuts 18/11/2008

if any(find(~isfinite(x+y)))
  error('no nans allowed in data')
end

aest = aest(:);


%%% normalize variables as appropriate
xoff = x(1);
xnorm = x - xoff;
xnormavg = mean(xnorm);

yavg = median(y);
ystd = std(y);
ynorm = y - yavg;


%%% set up fitting intervals with defaults, if necessary
if ~exist('int_lin','var') || isempty(int_lin)
  int_lin = [100 xnorm(end)];
end

if ~exist('int_exp','var') || isempty(int_exp)
  int_exp = [0 100];
end

if xnorm(end) < int_lin(1);
  warning(['x is not longer than 100 days - fitting may be ill ' ...
           'defined'])
  disp('using half-way value to set int_lin and int_exp')
  int_lin = xnorm([ round(length(xnorm)/2) end ]);
  int_exp = xnorm([ 0 round(length(xnorm)/2)]);
end


%%% remove the quadratic trend over int_lin (latter portion of data)
ii = find( xnorm >= int_lin(1) & xnorm <= int_lin(2) );

afit_quad1 = polyfit(xnorm(ii)-xnormavg,ynorm(ii),2);
yfit_quad1 = polyval(afit_quad1,xnorm-xnormavg);

ynorm2 = ynorm - yfit_quad1; % remove initial quadratic trend


%%% find int_exp region (initial portion of data) and fit exponential
ii = find( xnorm >= int_exp(1) & xnorm <= int_exp(2) );

h_expfun = inline('A(1) * exp( x/A(2) ) + A(3)','A','x');
afit_ex = nlinfit(xnorm(ii),ynorm2(ii),h_expfun,[aest([1:2])' 0]);
yfit_ex = h_expfun(afit_ex,xnorm);
yfit_ex2 = yfit_ex + yfit_quad1;

ynorm3 = ynorm - yfit_ex;  % remove exponential


%%% calculate final quadratic fit
afit_quad2 = polyfit(xnorm-xnormavg,ynorm3,2);
yfit_quad2 = polyval(afit_quad2,xnorm-xnormavg);

%keyboard


ynorm4 = ynorm - yfit_quad2;  % remove second quadratic fit


%%% calculate final exponential fit
afit_ex3 = nlinfit(xnorm(ii),ynorm4(ii),h_expfun,afit_ex);
yfit_ex3 = h_expfun([afit_ex3(1:2) 0],xnorm); % discard offset
yfit_ex4 = yfit_ex3 + yfit_quad2;


ynorm5 = ynorm - yfit_ex3;


%%% calculate third quadratic fit
afit_quad3 = polyfit(xnorm,ynorm5,2);
yfit_quad3 = polyval(afit_quad3,xnorm);

yfit_ex5 = yfit_ex3 + yfit_quad3;


%%% combine for final coefficients and final fit
yfit = yfit_ex5 + yavg;
yfit_quad = yfit_quad3 + yavg;


%%%% make sure mean offset is 0
%yavg2 = mean(y-yfit);
%yfit = yfit + yavg2;
%yfit_lin = yfit_quad + yavg2;


afit(1:2) = afit_ex3(1:2);
afit(3)   = afit_quad3(2);
afit(4)   = interp1(xnorm,yfit_quad,mean(xnorm));
afit(5)   = afit_quad3(1);

afit = afit(:);



if 1==0 % plot results

  if 1==0
    figure
    plot(xnorm,ynorm,'k-',xnorm,yfit_quad1,'b',...
         int_lin([1 1 1 2 2]),[ylim nan ylim],'.b:')

    plot(xnorm,ynorm2,'k-',xnorm,yfit_ex,'r',...
         int_exp([1 1 1 2 2]),[ylim nan ylim],'r:')
    hold on
    plot(xnorm,ynorm3,'b-')

    clf
    plot(xnorm,ynorm,'k-',xnorm,yfit_ex2,'r')

    plot(xnorm,ynorm3,'k-',xnorm,yfit_quad2,'r')

    plot(xnorm,ynorm4,'k-',xnorm,yfit_ex3,'r',...
         int_exp([1 1 1 2 2]),[ylim nan ylim],'r:')
    hold on
    plot(xnorm,ynorm4 - yfit_ex4,'b-')

    clf
    plot(xnorm,y,'k-',xnorm,yfit,'r',xnorm,yfit_quad,'b-')
    hold on
    plot(int_lin([1 1 1 2 2]),[ylim nan ylim],'or:',...
         int_exp([1 1 1 2 2]),[ylim nan ylim],'.b:')
  end

  figure(94), clf
  plot(xnorm,y,'k-',xnorm,yfit,'r',xnorm,yfit_quad,'b-')
  hold on
  plot(int_lin([1 1 1 2 2]),[ylim nan ylim],'.b:',...
       int_exp([1 1 1 2 2]),[ylim nan ylim],'or:')

end