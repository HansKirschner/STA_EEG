function [i,v,d] = closeval(x,y,varargin)
%returns the indices, closest values, and distance in x for every entry in y.
%uses min and thus returns the first closest value found if two have the
%same distance.
%m = 2 removes the entry from the possible set per loop iteration.
if isempty(varargin)
    m=1;
else
    m=varargin;
end

v = nan(length(y),1); i = v; d = v;
for c = 1 : length(y)
    [d(c),i(c)] = min(abs(x-y(c)));
    v(c) = x(i(c));
    if m == 2
        x(i(c)) = nan;
    end
end