function s=format_pval(p)
% just returns a p value formatted as a string if exp output in matlab is activated. 
%%
s = num2str(p);
pe = strfind(s,'e');
if ~isempty(pe)
    s2 = num2str(str2num(s(pe+2:end)));
    s = [s(1:4) 'x10^{-' s2 '}'];
end
return