function [m,mi]=meannan(data);
% function [m,mi]=meannan(data);
%
% calculates mean-values similar to built-in 'mean' 
% but takes NaN's into account
% by removing them and calculating the mean of the rest
%
% input  :	data		: data array
%
% output :	m		: mean values
%		mi		: mean of indices of the 'averaged' values
%
% version 1.1.0		last change 04.09.1995

% Gerd Krahmann, IfM Kiel, Mar 1993
% added mi-output	G.Krahmann, Sep 1995

b=isnan(data);
a=sum(b(:));
[s1,s2]=size(data);
indices=[1:s1]'*ones(1,s2);

if (a>0)
  if ( s1==1 ) + ( s2==1 )
    data=data(~b);
    indices=indices(~b);
    if isempty(data)
      m=nan;
      mi=nan;
    else
      m=mean(data);
      mi=mean(indices);
    end
  else
    m=zeros(1,s2);
    for i=1:s2
      [m1,m2]=mean_nan(data(:,i));
      m(i)=m1;
      mi(i)=m2;
    end
  end
else
  m=mean(data);
  mi=mean(indices);
end
