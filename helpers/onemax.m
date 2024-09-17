function x = onemax(x,W)
%finds the only one maximum (W = 1), minimum (W = 2), or extreme (W = 3)
%value in an n-dimensionaly array x.
%%
x = squeeze(x); %remove singletons
dims = size(x);
for c = 1 : numel(dims)
    if W == 1
        x = squeeze(max(x));
    elseif W == 2
        x = squeeze(min(x));
    elseif W == 3
        x = squeeze(max(abs(x)));
    else
        display('Unknown input for W'); return
    end
end

