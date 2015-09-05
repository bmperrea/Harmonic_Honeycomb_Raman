function [H,X,E]=histo(Y,bins)
% [H,X]=HISTO(Y,N)
% Calculates histogram with fixed binwidth N (default N=1). The used 
% method is very fast, exploiting advantages of sorted data
% HISTO(Y,N) 
%	plots the histogram bar(X,H)
% [H,X,E]=HISTO(Y,N) 
%	calculates also the Entropy E of Y with respect of binwidth N

%	Version 2.10
%	last revision 07.07.1997
%	Copyright (c) 1997 by Alois Schloegl
%	e-mail: a.schloegl@ieee.org	

% ROUND and SORT Data
Y=sort(Y); 

L=length(Y);
X=Y(1):Y(L);
H=zeros(size(X));

% transform range of Y to 1..max(Y)
%Y=Y+1-min(Y);

% find k where Y(k) is different to Y(k-1)
k=find(diff(Y)>0);
H(Y(k))=diff([k; length(Y)+1]);

%plots it
if nargout == 0
    	bar(X,H);
elseif nargout==3
% calculates ENTROPY = sum(p(i)*log(p(i))) with p(i)=H/M;
%p=H(Y(k))/L;
%E=-sum(p.*log(p))/log(2);
end;