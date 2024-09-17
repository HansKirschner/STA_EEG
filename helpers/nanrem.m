function [x,i]=ci_isfield(S,f)
% Case insensitive isfield function.
%%
value = [];
names   = fieldnames(S);
isField = strcmpi(f,names);  

if any(isField)
  value = S.(names{isField});
else
  
end
return