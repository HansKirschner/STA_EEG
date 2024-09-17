function [ x ] = make_wide( x )
%just makes sure a variable is long or wide
if size(x,1) > 1
    x=x';
end
return




















