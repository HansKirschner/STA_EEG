function [OD,AddInf] = STA_permutation(indata,s) 

if isfield(s.dimensions)
    dims = s.dimensions;
else
    dims = 2;
end

indata = rand(20,50); %data needs to be: subjects / trials, time, electrodes
n_iter = 1000;
n_factors = size(indata,1); 
n_permut = factorial(n_factors)/ (factorial(n_factors/2) * factorial(n_factors/2));

PermMatrix = randi([0 1],n_factors,n_iter).*2-1;

%extend matrices to be identical
%%
clc
% Matrix Vatiant uses more memory and is not faster than loop!
% clear A B C
% tic
% A = repmat(indata,1,1,n_iter);
% B = permute(repmat(PermMatrix,1,1,50), [1 3 2]);
% C = A.*B;
% toc

clear C
tic
C = nan([size(indata) n_iter]);
for c = 1 : n_iter
    C(:,:,c) = indata.*repmat(PermMatrix(:,c),1,50);
end
toc

%%
n_perm = 10000;
clear Longest_time Sum_T Mean_T d2
d1 = d(:,126:376);

for np = 1 : n_perm
    switcher = (2*round(rand(31,1))-1)';
    for pbn = 1:size(d1,1)
        d2(pbn,:) = d1(pbn,:) * switcher(pbn);
    end
    [h,~,~,stats]=ttest(d2);
    %plot([stats.tstat])
    start_c_p = 0; time_start_p = []; e_p_count=0;
    start_c_n = 0; time_start_n = []; e_n_count=0;
    time_of_effect_p=[];
    time_of_effect_n=[];
    for c = 1 : size(d1,2)
        if stats.tstat(c) > Tcrit & ~start_c_p %positive effect
            start_c_p = 1;
            time_start_p = c;
        end;
        if (stats.tstat(c) < Tcrit & start_c_p == 1) | (start_c_p == 1 & c == size(d1,2)) % No longer positive effect
            e_p_count = e_p_count + 1;
            start_c_p = 0;
            time_of_effect_p(e_p_count,1:2) = [time_start_p c];
            time_of_effect_p(e_p_count,3) = c-time_start_p;
            time_of_effect_p(e_p_count,4) = sum([stats.tstat(time_start_p:c)]);
            time_of_effect_p(e_p_count,5) = mean([stats.tstat(time_start_p:c)]);
        end;

         if stats.tstat(c) < -Tcrit & ~start_c_n %positive effect
            start_c_n = 1;
            time_start_n = c;
        end;
        if (stats.tstat(c) > -Tcrit & start_c_n == 1) | (start_c_n == 1 & c == size(d1,2)) % No longer positive effect
            e_n_count = e_n_count + 1;
            start_c_n = 0;
            time_of_effect_n(e_n_count,1:2) = [time_start_n c];
            time_of_effect_n(e_n_count,3) = c-time_start_n;
            time_of_effect_n(e_n_count,4) = abs(sum([stats.tstat(time_start_n:c)]));
            time_of_effect_n(e_n_count,5) = abs(mean([stats.tstat(time_start_n:c)]));
        end;
    end;
    Combi = [time_of_effect_p; time_of_effect_n ];
    if isempty(time_of_effect_n) & isempty(time_of_effect_p)
        Longest_time(np) = 0;
        Sum_T(np) = 0;
        Mean_T(np) = 0;
    else
        Longest_time(np) = max(Combi(:,3));
        Sum_T(np) = max(Combi(:,4));
        Mean_T(np) = max(Combi(:,5));
    end;
    if size(Combi,1) > 1
        Longest_time2(np) = max(Combi(2:end,3));
        Sum_T2(np) = max(Combi(2:end,4));
        Mean_T2(np) = max(Combi(2:end,5));
    else
        Longest_time2(np) = 0;
        Sum_T2(np) = 0;
        Mean_T2(np) = 0;
    end;
end;
%%
%get positions of ERN and Pe cluster compared to random samples
ERN_perc_length = length(Longest_time(Longest_time>ERN(3)))/5000
ERN_perc_sum    = length(Sum_T(Sum_T>ERN(4)))/5000
ERN_perc_mean   = length(Mean_T(Mean_T>ERN(5)))/5000

Pe_perc_length = length(Longest_time(Longest_time>Pe(3)))/5000
Pe_perc_sum    = length(Sum_T(Sum_T>Pe(4)))/5000
Pe_perc_mean   = length(Mean_T(Mean_T>Pe(5)))/5000


%Combined likelihood of observing two effects of at least equal size?
sum(ismember(find(Sum_T>ERN(4)), find(Sum_T2 > Pe(4))))