%%%ENTER HERE CRITERIA
function [  ] = STA_2nd_Level_Regression(infold, SLevelT, settings)
%%
%'infold'   = folder with first level output. Will generate a new folder
%               within here, where all second level outputs are stored.
%'slevelT'  = a table with fields corresponding to every first level
%               dataset.
%settings   with fields:
%   ModelName   =   name for the model (needs to be unique!)
%   exclude     =   vector of subjects to exclude (alphabetical order assumed, e.g. [4 5 6]
%   UseValues   =   'b' or 't' will use first level t or b values, 'EEG'
%                   will use returned EEG per bin from within subjects
%                   analysis (default: uses first returned value!)
%   EEGbin      =   bin from returned EEG data (only relevant if UseValues = 'EEG')
%   UseReg      =   vector of regressors to include (empty: all)

%%
ModelName   = {};   %Name for the second level model
exclude     = [];   %These participants will be excluded from the analysis
OutString   = '';   %A string to display information about the analysis
UseValues   = 'b';  %Default: Uses first level b values
UseReg      = [];   %Default: All regressors will be included for 2nd level models
RegLables   = {};   %lables for regression values (important for categorial regressors only, order low to high)
IncElect    = {};   %Electrodes to include (default: all, usage: {'FCz' 'Cz'})
RunParallel = 0;    %Default: Do not use parallel pool
EEGbin      = 1;    %Default, if set to use EEG, use first returned bin
outfold     = '';
Fdisp       = 1;
ERPimage    = 0;    %Default: No ERP image
ERPimageSm  = 0; 

if isfield(settings, 'ModelName');  ModelName       = settings.ModelName;       end
if isfield(settings, 'exclude');    exclude         = settings.exclude;         end
if isfield(settings, 'UseValues');  UseValues       = settings.UseValues;       end
if isfield(settings, 'UseReg');     UseReg          = settings.UseReg;          end
if isfield(settings, 'RegLables');  RegLables       = settings.RegLables;       end
if isfield(settings, 'IncElect');   IncElect        = settings.IncElect;        end
if isfield(settings, 'RunParallel');RunParallel     = settings.RunParallel;     end
if isfield(settings, 'EEGbin');     EEGbin          = settings.EEGbin;          end
if isfield(settings, 'outfold');    outfold         = settings.outfold;         end
if isfield(settings, 'ERPimage');   ERPimage        = settings.ERPimage;        end
if isfield(settings, 'ERPimageSm'); ERPimageSm      = settings.ERPimageSm;      end

if Fdisp > 0; fprintf('\n\n**********Running EEG Regression*********\n'); end

%Test settings and override defaults
if isempty(ModelName)
    disp(['Please specify a name for the second level model.']); return;
end
if ~strcmp(infold(end),'/')
    infold=[infold '/'];
end
if ~isempty(outfold) && ~strcmp(outfold(end),'/')
    outfold=[outfold '/'];
end

%Get files
F = dir([infold '*.mat']);
if isempty(F)
    disp(['No files found...']); return;
end
%
VPN = [1:length(F)];
if size(SLevelT,1) ~= length(VPN)
    disp('Designmatrix needs to be of the same length as the number of input datasets.');return
end

if ~isempty(exclude)
    OutString = [OutString 'Excluding ' num2str(length(exclude)) ' participants from the analysis.\n'];
    VPN = setdiff(VPN, exclude);
    SLevelT = SLevelT(VPN,:);
end
Ninclu = size(VPN,2);

save_path = [infold '2ndlevel/' outfold  ModelName '/'];
if ~exist(save_path)
    mkdir(save_path);
    disp(['Creating new folder for model ' outfold ModelName])
end

if RunParallel
  parforArg = Inf;
else
  parforArg = 0;
end


%load first subjects and termine fields of output
VP              = load([infold F(1).name]);
FN1             = fieldnames(VP);
FN2             = fieldnames(VP.(FN1{1}));
FN3             = fieldnames(SLevelT);
Info            = VP.(FN1{1}).(FN2{2});
Data            = VP.(FN1{1}).(FN2{1});
Electrodes      = Info.Output_Labels;
Elocs           = Info.Output_Electrodes;
dims            = size(Data.t_values);
RegressionTime  = Info.Return_Timewindow_ms;
MulM            = 1; %Dummy: can be changed by user of more than one model is there

%channel location information
try
    A = settings.chanlocs;
catch
    A = load('data/ChanLocs.mat');
end

AF              = fieldnames(A);
chanlocs        = A.(AF{1});
ChanLoc         = chanlocs(Elocs);

%adjust inputs
try Info.IncluReg
catch
    Info.IncluReg = 1 : length(Info.RegNames);
end
if isempty(UseReg)
    UseReg = 1 : length(Info.RegNames);
end
    
%Check which regressors should be included
RegPresent  = Info.IncluReg;
UseReg      = intersect(Info.IncluReg, UseReg);

%check which electrodes should be included
if ~isempty(IncElect)
    IncElect = structfind(ChanLoc,'labels',IncElect);
else
    IncElect = structfind(ChanLoc, 'labels', {chanlocs(Elocs).labels});    
end

%pre allocate arrays
dims(1)     = length(IncElect); %reduce data size to actual number of included electrodes
dims(2)     = length(UseReg);   %reduce data size to actual number of included regressors
try TF = Info.TF;
catch
    try TF = Data.TF;
    catch
        TF = 'Time domain';
    end
end

if isstruct(TF)
    first_level  = nan(Ninclu, dims(1), dims(2), dims(3), dims(4)); %5 dims = VP, electrode, regressors, time, frequency
    OutString = [OutString 'Uses time-frequency data.\n'];
else
    first_level = nan(Ninclu, dims(1), dims(2), dims(3)); %4 dims = VP, electrode, regressors, time
    OutString = [OutString 'Uses time-domain data.\n'];
end

%check which values should be used
switch UseValues
    case 'b'
        OutString = [OutString ' Using first level b values.\n'];
    case 't'
        OutString = [OutString ' Using first level t values.\n'];
    case 'EEG'
        OutString = [OutString ' Using first level EEG from bin ' num2str(EEGbin) '.\n'];
        for cEEGBin = 1 : length(UseReg); OutString = [OutString '\tThis corresponds to regressor ' Info.RegNames{cEEGBin} ' field: ' Info.RegLables{UseReg(cEEGBin),EEGbin} '\n'];end
    otherwise
        disp(['Unknown value for field settings. UseValues: ' settings.UseValues]); return
end

%%
%Generate display output
clc
disp(['Found ' num2str(Ninclu) ' files in folder ' infold '.'])
disp(['Regression model ' strrep(FN2{1},'_',' ') ' includes ' num2str(length(IncElect)) ' electrodes and ' num2str(min(Info.Return_Timewindow_ms)) ' until ' num2str(max(Info.Return_Timewindow_ms)) 'ms.'])
fprintf([OutString '\n'])
fprintf(['\nNumber of orignial regressors is ' num2str(length(Info.RegNames)) ':\n'])
for c = 1 : length(Info.RegNames)
    if ismember(c, UseReg)
        fprintf(['\tReg no ' num2str(c) ' (included): ' Info.RegNames{c} '\n'])
    else
        fprintf(['\tReg no ' num2str(c) ' (excluded): ' Info.RegNames{c} '\n'])
    end
end

fprintf(['\n\nDesign matrix for second level contains fields: \n'])
for c = 1 : length(FN3)-3
    Uval = unique(SLevelT.(FN3{c}));
    fprintf(['\tReg no ' num2str(c) ': ' FN3{c} ' (' num2str(length(Uval)) ' unique values)\n'])
    if isempty(RegLables) %no lables for values provided, generating new ones
        if length(Uval)>2
            NewRegLables(c,1) = {'low'};
            NewRegLables(c,2) = {'high'};
        elseif length(Uval)>1
            NewRegLables(c,1) = {'numeric'};
            NewRegLables(c,2) = {'numeric'};
        else
            disp('Warning: Designmatrix contains only one value... Constant will be automatically added! Aborting.'); return;
        end
    else
        NewRegLables(c,1) = {RegLables{c,1}};
        NewRegLables(c,2) = {RegLables{c,2}};
    end
    fprintf(['\t\tValues: ' NewRegLables{c,1} ' and ' NewRegLables{c,2} '.\n'])
end
%%
fprintf('Loading data from: ')
for p = 1 : Ninclu
    %Get fitting INFO file and display it
    VPname = [F(VPN(p)).name(1:end-4) '  '];
    if p > 1 && p~=Ninclu
        fprintf(2,[RemoveString '%s'], VPname);
    elseif p == Ninclu
        fprintf(2,[RemoveString '%s'], VPname);
        fprintf('\n\n')
    else
        fprintf(2,VPname);
    end
    RemoveString = [VPname VPname]; RemoveString(1:2:length(RemoveString))='\';RemoveString(2:2:length(RemoveString))='b';
    %load actual data
    VP = load([infold F(VPN(p)).name]);
    if isstruct(TF)
        if strcmp(UseValues,'t')
            first_level(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).t_values(IncElect,UseReg,:,:); %5 dims = VP, electrode, regressors, time, frequency
        elseif strcmp(UseValues,'b')
            first_level(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).b_values(IncElect,UseReg,:,:); %5 dims = VP, electrode, regressors, time, frequency
        else
            disp(['Unrecognized value for parameter UseValues: ' UseValues]);return;
        end
    else
        if strcmp(UseValues,'t')
            first_level(p,:,:,:) = VP.(FN1{1}).(FN2{MulM}).t_values(IncElect,UseReg,:,:);
        elseif strcmp(UseValues,'b')
            first_level(p,:,:,:) = VP.(FN1{1}).(FN2{MulM}).b_values(IncElect,UseReg,:,:);
        elseif strcmp(UseValues,'EEG')
            %%
            DummyEEG = permute(VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UseReg,EEGbin,IncElect,:), [3 1 2 4]); %bring in regression order = Elec, Regressor, Time
            s = size(DummyEEG);
            DummyEEG = reshape(DummyEEG,[s(1) s(2) s(4) 1]); %remove dimension for EEG bin
            first_level(p,:,:,:) = DummyEEG;
            %%
        else
            disp(['Unrecognized value for parameter UseValues: ' UseValues]);return;
        end
    end
end
%%
%time-frequency regression
if isstruct(TF)
    TF.SecondLevel = 1;
    % parfor (c = 1:length(UseReg), parforArg)
    for c = 1:length(UseReg)
        %loop through regressors
        Reg2=struct();
        Reg2.chanlocs = ChanLoc(IncElect);
        Reg2.times = Info.Return_Timewindow_ms;
        %Add fields to make this EEG structure like
        DummyEEG = permute(first_level(:,:,c,:,:), [2 4 1 3 5]); %Bring into EEGlab / regression order (Electrode, Time, Epoch/Subject, frequency)
        s = size(DummyEEG);
        Reg2.data = reshape(DummyEEG,[s(1) s(2) s(3) s(5) 1]); %removes third dimension
        %Run Regression
        [ Odata, add_info ] = STA_Fast_Regress(Reg2, [1:length(IncElect)], [Reg2.times(1) Reg2.times(end)], 'bin_size', 0, ...
             'stepsize',1, 'TF', [TF], 'model 1', {{1:length(VPN)} SLevelT {'EEG'} {'Level2Reg1'}},'binEEG', 3, 'retall', 3, ...
             'PredNames', {FN3{1:end-1}}, 'PredLabels', RegLables, 'normaliseon', 0, 'Regmode', 'ols');
        OD(c) = {Odata}; AI(c)={add_info};
    end

    for c = 1 : length(UseReg)
        clear OutputFirstLevel
        OutputFirstLevel.(Info.RegNames{UseReg(c)}) = OD{c};
        OutputFirstLevel.Info = AI{c};
        OutputFirstLevel.Info.ModelName = ModelName;
        OutputFirstLevel.Info.Output_Electrodes = IncElect;
        OutputFirstLevel.Info.UseValues2ndLevel = UseValues;
        OutputFirstLevel.Info.EEGbin2ndLevel = EEGbin;
        save([save_path ModelName ' ' Info.RegNames{UseReg(c)}], 'OutputFirstLevel')
    end
end
%%
%time domain regression
if ~isstruct(TF)
%     parfor (c = 1:length(UseReg), parforArg)
    for c = 1:length(UseReg)
        %loop through regressors
        Reg2=struct();
        Reg2.chanlocs = ChanLoc(IncElect);
        Reg2.times = Info.Return_Timewindow_ms;
        %Add fields to make this EEG structure like
        Reg2.data = permute(first_level(:,:,c,:), [2 4 1 3]); %Bring into EEGlab order (Electrode, Time, Epoch/Subject)
        if ~isfield(Reg2,'srate')
            Reg2.srate = 500;
        end
        %Run Regression
        [ Odata, add_info ] = STA_Fast_Regress(Reg2, [1:length(IncElect)], [Reg2.times(1) Reg2.times(end)], 'bin_size', 0, ...
             'stepsize',1, 'TF', [], 'model 1', {{1:length(VPN)} SLevelT {'EEG'} {'Level2Reg1'}},'binEEG', 3, 'retall', 3, ...
             'PredNames', {FN3{1:end-1}}, 'PredLabels', RegLables, 'normaliseon', 0, 'Regmode', 'ols','ERPimage',ERPimage,'ERPimageSmooth',ERPimageSm);
        OD(c) = {Odata}; AI(c)={add_info};
        
%     end
% 
%     for c = 1 : length(UseReg)
        clear OutputFirstLevel
        OutputFirstLevel.(Info.RegNames{UseReg(c)}) = OD{c};
        OutputFirstLevel.Info = AI{c};
        OutputFirstLevel.Info.ModelName = ModelName;
        OutputFirstLevel.Info.Output_Electrodes = IncElect;
        OutputFirstLevel.Info.UseValues2ndLevel = UseValues;
        OutputFirstLevel.Info.EEGbin2ndLevel = EEGbin;
        if strcmp(UseValues,'EEG')
            OutputFirstLevel.Info.EEG_Factor_Level = Info.RegLables{UseReg(c),EEGbin};
        else
            OutputFirstLevel.Info.EEG_Factor_Level = '';
        end
        
        save([save_path ModelName ' ' Info.RegNames{UseReg(c)}], 'OutputFirstLevel')
        
        
        
    end
end


















