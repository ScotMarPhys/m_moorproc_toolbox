function [out]=filtnan(b,a,data);
% function [out]=filtnan(b,a,data);
%
% similar to 'filtfilt.m' but is able to handle NaN's
% NaN's inside other data are linearly interpolated
% NaN's at the end are replaced by last value
%
% input  :	- b	: order of filtering
%		- a	: cut-off frequency 
%		- data  : unfiltered data
%
% output :	- out	: filtered data
%
% version 1.0.1		last change 15.03.1996

% G.Krahmann, IfM Kiel, Sep 1994

bad=find(~isnan(data));
if ~isempty(bad)
  s=size(data);
  if (s(2)==1)
    data=data';
    turn=1;
  else
    turn=0;
  end
  l=length(data);
  good=find( ~isnan(data) );
  if (good(1)~=1)
    data(1:good(1)-1)=data(good(1))*ones(1,good(1)-1);
  end
  good=find( ~isnan(data) );
  if (good(length(good))~=l)
    data(good(length(good))+1:l)=...
	data(good(length(good)))*ones(1,l-good(length(good)));
  end
  bad=find( isnan(data) );
  if ~isempty(bad)
%    size( ([1:l]<bad(1))  )
    for n=1:length(bad)
      i1=max( find( (~isnan(data)) & ([1:l]<bad(n)) ) );
      i2=min( find( (~isnan(data)) & ([1:l]>bad(n)) ) );
      data(bad(n))=data(i1)+(data(i2)-data(i1))/(i2-i1)*(bad(n)-i1);
    end
  end
  if (turn==1)
    data=data';
  end
end
out=filtfilt(b,a,data);
