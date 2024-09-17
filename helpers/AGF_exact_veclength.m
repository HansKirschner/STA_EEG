function [ Out ] = AGF_exact_veclength( S, E, L)
%Function returns a vector from S to E with the exact length L.
%%
%Determine correct step size to get an array of exactly 100
%L = 100;
%E = 200;
StepSize = E / L; %start is close to usual stepsize, but this can vary (maybe there is a better solution??)
stepcor = 0;  
while ~stepcor
    Vec = size([S : StepSize : E],2);
    if Vec == L
        stepcor = 1;
    elseif Vec < L
        StepSize = StepSize * 0.99;
    else
        StepSize = StepSize * 1.01;
    end;
end;
Out = [1 : StepSize : E];
