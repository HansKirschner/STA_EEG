function [rej_trials] = STA_adaptive_artifact(EEG, TI, S)
%This function iteratively uses pop_jointprob to remove a desired target
%number of trials. Trials are defined by the trial index (TI), such that no
%more than a certain number of trials can be set to be removed. 
%AGF, 2015
%%
%adaptive artifact rejection
%Default settings:
max_p   = 10;        % maximum percentage of trials allowed to be removed by artifact rejection
min_n   = 3;         % at least n trials should be rejected -> can be set, but not necessarily so
StartSD = 5;         % starting position for optimization
EpLimit = [];        % limits to re-epoch of wanted (unit = seconds!)
StepSize = 0.1;      % step size for increasing or decreasing threshold

%overwrite with settings (S)
if isfield(S, 'max_p');     max_p = S.max_p;     end
if isfield(S, 'min_n');     min_n = S.min_n;     end
if isfield(S, 'StartSD');   StartSD = S.StartSD; end
if isfield(S, 'EpLimit');   EpLimit = S.EpLimit; end

if islogical(TI) %convert logical into index
    TI = find(TI);
end

EEG = AGF_epoch_by_number(EEG, TI, [], EpLimit); %reduce trials

disp(['Removing at least ' num2str(min_n) ' trials and not more than ' num2str(max_p) '% of trials (n = ' num2str(length(TI)) ').'])

%%
Fertig      = 0;        % while loop index
UseThresh   = StartSD;  % Initial SD 

while ~Fertig
    [T1, ~, ~, nrej] = pop_jointprob_silent(EEG,1,1:size(EEG.data,1),UseThresh,UseThresh,0,0);
    if nrej/size(EEG.epoch,2) > max_p/100        %too many removed
        UseThresh = UseThresh + (UseThresh*StepSize);    %Increase threshold by StepSize
        disp(['Increasing threshold by ' num2str(StepSize*100) '%. Current threshold is: ' num2str(UseThresh)])
    elseif nrej < min_n                         %not enough have been removed
        UseThresh = UseThresh - (UseThresh*StepSize);    %Decrease threshold by StepSize
        disp(['Lowering threshold by ' num2str(StepSize*100) '%. Current threshold is: ' num2str(UseThresh)])
    else
        Fertig=1;
    end;
end;
disp([num2str(nrej) ' trials removed.'])
rej_trials = TI(find(T1.reject.rejjp==1)); %trials to reject in index of initial trials
return;
