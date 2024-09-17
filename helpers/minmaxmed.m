function [v,t] = minmaxmed(d,m)
% Returns min, max, mean or median of d depending on 'm'.
%%
v=[];t=[];
if strcmpi(m, 'mean')
    v = mean(d);
elseif strcmpi(m, 'median')
    v = median(d);
elseif strcmpi(m, 'min')
    [v,t] = min(d);
elseif strcmpi(m, 'max')
    [v,t] = max(d);
end

return