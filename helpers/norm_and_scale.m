function x = norm_and_scale(x,r1,r2)
% Normalise values of an array to be between r1 and r2
% original sign of the array values is maintained.
x = (x-min(x))*(r2-r1)/(max(x)-min(x)) + r1;
return;
