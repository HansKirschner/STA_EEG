function v = ci_isfield(S,f)
% Case insensitive isfield function.
%%
v       = [];
names   = fieldnames(S);
isField = strcmpi(f,names);  

if any(isField)
    v = S.(names{isField});
end
 return