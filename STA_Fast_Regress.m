function [ Odata, add_info ] = STA_Fast_Regress( indata, electrodes, time_w, varargin)
%%
%Function that regresses EEG data.
%AGF HK 2022
%ToDo:
% - overwork included single-trial baseline correction (currently deactivated)

%Added option to define a range of wavelet cycles (e.g. [3 13]), so that the width changes
% as a function of frequency
%Added information for normalization of ERPs when plotting, added option to
%   exclude regressors from output (will still be included in model)
%Added single-trial baseline correction for TF (05/18/2018)
%Added length comparison for labels and names (6/10/2015)
%Fixed when removing NaNs from design matrix (10/06/2015)
%Fixed bugs when using EEG as predictor with new lables and switched dimensions of Y (10/05/2015)
%--
%Returns more information about regressors including numbers per items for SE calculation on second level etc. (09/23/2015)
%Now also has a string containing warnings for easier bugfixing.
%--
%Added compression of ERPimage data (default: on - set ImageCompression = 0 to turn off) by rounding and converting to single. Warning: This may be imprecise and is thought for plotting only! (09/14/2015)
%--
%Now uses ols script (from Tim Behrens) to run regression. Also supports automatic design matrix changes (08/09/2015)
%--
%Include sorted ERP image plots and changed loop to only calculate binned EEG once per electrode (08/05/2015)
%--
%Uses less memory for TF analysis (05/12/2015)
%--
%Corrected R2 calculation dimension bug (05/18/2015)
%INPUT:
%'indata'       either EEG structure or path to eeg file that will be
%               openend.
%'electrodes'   Electrodes numbers to be included.
%'time_w'       time window (ms) that should be used for regression analyses. 
%
%OPTIONAL INPUT:
%
%'bin_size'     Smoothing over this amount of datapoints before and after the current datapoint.
%-
%'stepsize'     Performs the test in every n-th datapoint (default = length(time_w) -> only one test.
%-
%'RegMode'      'Robust' uses matlab robustfit (slower), 'ols' uses faster 'regress', 'pinv' uses pseudo inverse (fastest)
%-
%'NormaliseOn'  If (1) will normalise predictor matrix before regression. 
%-
%'TF'           Structure with settings for TF analysis with fields:
%                   TF.frequencies  = frequency range
%                   TF.stepnumber   = number of steps within that range
%                   TF.cyclenumber  = number of cycles (will scale with the frequency (FWHM of gaussian) 
%                                     can also be a range (e.g. [3 13]), so that the width changes
%                                     as a function of frequency
%                   TF.space        = 'linear' or 'log'
%                   TF.bands        = name for bands over which should be
%                                     averaged before(!) regression {'alpha' 'beta'}
%                   TF.bandfreq     = frequencies for those bands (e.g. {[8 12] [12 23]}), length must equal TF.bands field.
%                   TF.basetime     = baseline used for TF subtraction (most be included in the time-range!)
%                   TF.transform     = what happens with the data after TF transformation?
%                       'none'            - do absolutely nothing (or leave empty)
%                       'normalize'       - normalize power within each frequency
%                       'log'             - just take the log of the frequency power
%                   TF.basetype
%                       'subtract'        - subtracts each individual trials� baseline within each frequency
%                       'percent'         - subtracts each individual trials� baseline and divides by the mean baseline over all trials (determined prior to subtraction)

%-
%'base'         Specify certain periods per single trial as baseline (has
%               to be a 2D array of length (trial_ind) where column one =
%               first timepoint and column 2 = second timepoint.
%-
%'model'        Provide a model to perform robust regression in. Model has
%               to have the following structure = {{trialindex} {predictors} {Y} {name}}
%               trialindex = trials to be included relative to trials in the dataset (i.e., EEG.epoch).
%               For example if trials 1 4 7 are indexed, the regressor 1 0 1 would be relative to the index, 
%               not the initial trials 1 2 3. 
%               {name} has to be a string that names the model.
%               If {Y}='EEG' regression will use EEG as observation in the model.
%               If {Y}= any numeric regressor, EEG will be added as a predictor to the model (thus, EEG predicts another factor, not vice versa).
%-
%'PredNames'    Can be a cell with names for every regressor {'Reg1' 'Reg2'}.
%-
%'PredLables'   Can be a cell array of lables for categorial regressors (leave empty for parametric regressors!).
%               Note: The labels have to be in ascending order for their corresponding value. E.g. if 0 = 'male' and 1 = 'female' then
%               PredLables have to be {'male' 'female'}.
%-
%'binEEG'       Also returns regressors binned in high and low values (if
%               parametric) or per condition. 1 or 2 = median split, 3 = three
%               bins etc.
%               Output is field 'EEG_per_regressor' = 4D with Regressor, bin, electrode, time.
%               Field 'EEG_vlaues' stores the values per condition in each regressor (parametric = borders between bins).
%-
%'ERPimage'     Adds an ERP image to the output. Note: If this should later be used for averaging, set ERPimageDim to 
%               a defined value (e.g. 100). Images are then scaled to this value. 
%               'ERPimageDim' sets dimensions (default = []).
%               'ERPimageSmooth' sets smoothing (default = 20).
%               Output is stored as 4D field 'ERP_Image' with regressor | electrode | time | trials.
%               Info field 'add_info.ERPimage_Regressor_no' stores for which regressors an image has been calculated.
%-
%'NewEpoch'     Struct with fields for new epoching, which warps the output to this event:
%                 Time          = time in relation to current event (per trial) as a vector
%                 boundaries    = new epoch boundaries vector of times in ms, e.g., [-300 850]
%'RetAll'       0 = only returns b and t-values, 1 = returns pseudo R2 additionally, 2 = returns R2 for each Regressor, 3 = also returns p value time-course
%'IncluReg'     Disregard all other regressors (e.g., if IncluReg = [2 4]), will only return results for regressors 2 and 4 but accounts for the other factors nonetheless.
%'UseIca'       Cell array . If {1} will run regression on only this IC, if {1; [1 3 5]} will use activity of IC 1 for first analysis,
%               and [1 3 5] for next. 
%'IcaName'      Name you ICs (otherwise will use IC numbers). 
%'IcaCSpace'    1 = remain in IC space, 0 = use backprojection. If 1, correlates activity at this electrode with component activation. If this is negative, flips component activation. 
%'IcaElectrode' Electrode for backprofejction (if CPsace == 0), or for correlation to determine polarity (== 1). 
%REQUIRED FUNCTIONS
%   erpimage_noplot.m   For ERP images.
%   fconv_tg.m          For time fequency analyses.
%   ols.m               For faster ols regression.
%   normalise.m
%   nanrem.m

%Set defaults
%get datapoints from milliseconds
TimeInDP        = [find(indata.times==time_w(1)):find(indata.times==time_w(2))];
TimeInMS        = indata.times(TimeInDP);
TF              = [];              % Default: Work in time domain
Fmodel          = {};              % Default: No model used, no single datapoint regression.
nmodels         = 0;               % Default: No additional models are specified
Fdisp           = 1;               % Default: Display is activated
Fperc_base      = 0;               % Default: Do not perform transposition into % baseline activity
Fbase           = 0;               % Default: No baseline specified
FretEEG         = 0;               % Default: Do not return binned EEG per regressor
FStepSizeI      = NaN;
FbinsI          = 0;
NewEpoch        = [];
RetAll          = 0;               % Default: Do not return R2 statistics and p values
ERPimage        = 0;               % Default: No ERP image is returned
ERPimageDim     = 100;             % Default: ERP image has dimension of 100 (otherwise it is difficult to compare across participants!)
ERPimageSmooth  = 20;              % Default: Smoothing trials for ERP image is 20
RegMode         = 'robust';        % Default: Use robust regression.
NormaliseOn     = 1;               % Default: Normalises all predictors and replaces 0 with -1 values for dichotomous regressors.
ImageCompression= 1;               % Default: Compress ERP images.
Level           = 1;               % Default: First level analysis.
Warnings        = {};              % Empty Cell array to store warning messages.
PredNames       = {};              % Empty Cell array to store Names for predictors.
PredLables      = {};              % Empty Cell array to store values of predictors, i.e., their content.
ReturnInverted  = 0;               % Switch this to 1 if when using EEG as predictor, all other variables should be returned as well.
TF_Fdisp        = 2;
IncluReg        = [];
DownSample      = 1;               % Downsample data by this factor before returning (e.g., 2 = every second datapoint, 3 every third etc.)
FUseIca         = {};              % Default: Use raw data, not IC components.
FIcaName        = {};              % Default: No IC names.
FIcaElectrode   = 'Cz';            % Default: Electrode used in ICA backprojection
FIcaCSpace      = 0;               % Default: Use backprojection
if Fdisp > 0; fprintf('\n\n**********Running EEG Regression*********\n'); end
%%
for hide = 1:1 %Get Variables from varargin array
    nargs = nargin-3;
    if nargs > 1 
      if ~(round(nargs/2) == nargs/2)
        error('Odd number of input arguments??')
      end
    end
    for i = 1:2:length(varargin)
        Param = varargin{i};
        if ~isstr(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
        switch Param
            case 'perc_base'
                Fperc_base=lower(varargin{i+1});
                if isnumeric(Fperc_base)
                    if Fperc_base ~= 0 & Fperc_base ~= 1
                        disp(['Unrecognized value ' num2str(Fperc_base) ' for baseline input. has to be either 1 (active) or 0.']); return;
                    end
                else
                    disp(['Unrecognized value ' Fperc_base ' for baseline input. has to be either 1 (active) or 0.']); return;
                end
            case 'base'
                Fbase=lower(varargin{i+1});
                if size(Fbase,2)~=2
                    disp('Error: Baseline matrix has to be 2D with baselines over trials in columns...'); return;
                end
            case 'regmode'
                RegMode=varargin{i+1};
            case 'normaliseon'
                NormaliseOn=varargin{i+1};
            case 'bin_size'
                FbinsI=varargin{i+1};
            case 'bineeg'
                FretEEG=varargin{i+1};
            case 'stepsize'
                FStepSizeI=varargin{i+1};
            case 'disp'
                Fdisp=varargin{i+1};
            case 'tf'
                TF=varargin{i+1};
            case 'retall'
                RetAll=varargin{i+1};
            case 'newepoch'
                NewEpoch=varargin{i+1};
            case 'erpimage'
                ERPimage=varargin{i+1};
            case 'erpimagedim'
                ERPimageDim=varargin{i+1};
            case 'erpimagesmooth'
                ERPimageSmooth=varargin{i+1};
            case 'prednames'
                PredNames=varargin{i+1};
            case 'predlabels'
                PredLables=varargin{i+1};
            case 'inclureg'
                IncluReg=varargin{i+1};
            case 'downsample'
                DownSample=varargin{i+1};
            case 'useica'
                FUseIca=varargin{i+1};
            case 'icaname'
                FIcaName=varargin{i+1};
            case 'icaelectrode'
                FIcaElectrode=varargin{i+1};
            case 'icacspace'
                FIcaCSpace=varargin{i+1};
            case ['model ' num2str(nmodels+1)]
                nmodels=nmodels+1;
                Fmodel{nmodels}=varargin{i+1};
                if size(Fmodel{nmodels},2)~=4
                    disp('Models have to have 4 entrys: 1st = trials, 2nd = predictors, 3rd = observations,  4th = name. Leave 2nd empty {} to just use EEG.');return;
                end
            case 'analysislevel'
                Level=varargin{i+1};
                if Level == 1 && Fdisp > 0
                    disp('Running within subject regression.')
                elseif Level == 2 && Fdisp > 0
                    disp('Running second level regression.')
                else
                    disp('Unknown setting for input variable ''Level''...'); return;
                end
            otherwise
                disp(['Unknown argument ' varargin{i} '...'])
                pause
        end
    end 
    %Check input
    if strcmp(RegMode,'ols') 
        if Fdisp > 0
            disp('Using OLS regression.')
        end
    elseif strcmp(RegMode,'robust')
        if Fdisp > 0
            disp('Using robust regression.')
        end
    elseif strcmp(RegMode,'pinv')
        if Fdisp > 0
            disp('Using pseudo inverse solution.')
        end
    else
        disp(['Argument for RegMode ' RegMode ' not recognized.']); return
    end
    if FIcaCSpace~=1 && FIcaCSpace~=0;          disp(['Input argument ''IcaCSpace'' has to be logical.']); return; end
    if iscell(electrodes);                      electrodes  = structfind(indata.chanlocs,'labels',electrodes); end
    if isempty(electrodes) && isempty(FUseIca); error('Electrodes must be specified when not using ICA based regression.'); end
        
    
end

%%
%Set simple variables and pre allocate arrays
Fbins           = 1; %length(TimeInDP);                 % Default: Average over length of time-window
FStepSize       = 1; %length(TimeInDP);                 % Default: Average over length of time-window
Trials          = Fmodel{1}{1}{:}'; 
if ~isnan(FbinsI)
    Fbins=FbinsI;
end
if ~isnan(FStepSizeI)
    FStepSize=FStepSizeI;
end

disp(['Bin size is: ' num2str(Fbins) ' | Stepsize is : ' num2str(FStepSize)])
%Add some initial information to output
%add_info.Statistic_Timewindow_ms=[indata.times(TimeInDP(1)) indata.times(TimeInDP(end))];
%add_info.Statistic_times_ms = [indata.times(TimeInDP(1)):2:indata.times(TimeInDP(end))];
if isfield(indata,'srate')
    add_info.Bin_size_ms=Fbins*1000/indata.srate;
    add_info.Step_size_ms=FStepSize*1000/indata.srate;
else
    add_info.Bin_size_ms=Fbins;
    add_info.Step_size_ms=FStepSize;
end
add_info.Scale_2_Baseline_Percent=Fperc_base;
add_info.EEG_is_predictor=0;
%%
if isfield(indata, 'icaact') %reconstruct ica activity if necessary
    if isempty(indata.icaact) 
        indata.icaact = eeg_getica(indata);
    end
end
% keyboard
%% Replace data with IC activity
if ~isempty(FUseIca)
    n_comps     = length(FUseIca);  % number if indpendent components to be used
    electrodes  = 1 : n_comps;      % just replace electrodes with ICs 
    ICA_Electrode = structfind(indata.chanlocs,'labels',{FIcaElectrode});
    
    if isempty(ICA_Electrode)
        ICA_Electrode = 1;
    end
    
    for c = 1 : n_comps
        BPcomponent = FUseIca{c};
        if iscell(BPcomponent); BPcomponent=BPcomponent{:}; end
        %get correct labels
        try OutLable{c} = FIcaName{c};
        catch 
            if length(FUseIca{c}) == 1; OutLable{c} = ['IC_' num2str(FUseIca{c})];
            else
                ICNStr = 'ICs';
                for c2 = 1 : length(FUseIca{c})
                    ICNStr = [ICNStr '_' num2str(FUseIca{c}(c2))];
                end
                OutLable{c}=ICNStr;
            end
        end
        indata.chanlocs(c).labels = OutLable{c};
        %reduce data to ic activity
        
        if FIcaCSpace %use component activation
            for cc = 1 : length(BPcomponent)
                ICact(cc,:,:)       = squeeze(indata.icaact(BPcomponent(cc),:,:));  %extract IC activity for this component
                EEGact              = squeeze(indata.data(ICA_Electrode,:,:));      %and EEG activity
                r                   = corr([reshape(ICact(cc,:,:), [], 1) reshape(EEGact, [], 1) ]); %Get correlation between both
                ICact(cc,:,:)       = ICact(cc,:,:).*sign(r(2));                    %flipped IC activity 
                ICwinv(cc,:)        = indata.icawinv(:,BPcomponent(cc)).*sign(r(2)); %flipped topography
                if sign(r(2)) == 1
                    disp(['Using unflipped activity of component ' num2str(BPcomponent(cc))])
                else
                    disp(['Using flipped activity of component ' num2str(BPcomponent(cc))])
                end
                add_info.Corr_IC_EEG(cc) = r(2);
            end
            if cc > 1
                T2(c,:,:)           = squeeze(mean(ICact,1));
                add_info.ICA_Topo(c)= {squeeze(mean(ICwinv))};
            else
                T2(c,:,:)           = squeeze(ICact);
                add_info.ICA_Topo(c)= {squeeze(ICwinv)};
            end
        else %use backprojection
            disp(['Backprojection for: ' num2str(BPcomponent) ' at electrode: ' FIcaElectrode])
            T                   = pop_subcomp(indata, setdiff(1:size(indata.icawinv,2), BPcomponent), 0); 
            T2(c,:,:)           = T.data(ICA_Electrode,:,:);
            add_info.ICA_Topo(c)= {zscore(T.data(:,1,1))}; %get topography of the backprojected data
        end
        
    end
    %replace EEG data with pure backprojection of these components
    indata.data(1:n_comps,:,:)  = T2;
    indata.data(n_comps+1:end,:,:) = nan;
    clear T T2
    add_info.ICA_Used           = 1;
    add_info.ICs_Included       = FUseIca;
else
    add_info.ICA_Used           =0;
    add_info.ICs_Included       = [];
end
add_info.ICAspace           = FIcaCSpace;
add_info.ChannelLocations   = indata.chanlocs;
add_info.Output_Electrodes  = electrodes;

%%
for c = 1:length(electrodes) 
    if isfield(indata, 'chanlocs') && isempty(FUseIca)
        OutLable{c}=indata.chanlocs(electrodes(c)).labels;
    elseif isempty(FUseIca)
        OutLable{c}=['Electrode_' num2str(c)];
    end
end
add_info.Output_Labels=OutLable;  
    
if Fbase
    if size(Fbase,1)==1   
        FUbase=repmat(Fbase,size(indata.data,3),1);
    else
        FUbase=Fbase;
    end
    for trial = 1 : size(indata.data,3)
        if mod(FUbase(trial,1),2)
            b1=find(indata.times==FUbase(trial,1)-1);
        else
            b1=find(indata.times==FUbase(trial,1));
        end
        if mod(FUbase(trial,2),2)
            b2=find(indata.times==FUbase(trial,2)-1);
        else
            b2=find(indata.times==FUbase(trial,2));
        end
        %Get single trial baseline
        bsl=repmat(mean(indata.data(:,b1:b2,trial),2),1,size(indata.data,2));
        %subtract by baseline
        Data_base(:,:,trial)=indata.data(:,:,trial)-bsl;
    end
    add_info.Individual_Trial_Baseline=[mean(FUbase(:,1)) mean(FUbase(:,2))];
else
    add_info.Individual_Trial_Baseline=NaN;
end

%Perform security check: Is enough data at the borders of the time-window for the selected bin size?
if TimeInDP(1)-ceil(Fbins/2) <= 0 | TimeInDP(end)+floor(Fbins/2) > size(indata.data,2)
    disp('Warning! Timewindow + bin size to average is bigger than input data... stopping.')
    return;
end

%are labels and regressor lengths equal? 
if isempty(PredNames)
    disp('No names for predictors provided.')
elseif ~isempty(PredNames) && ~isempty(PredLables)
    
    if length(PredNames) ~= size(PredLables,1)
%        warning('Length of predictors and labels does not match, check for mistake?');pause(2);
    end
end

%Run the regression analysis
if isfield(indata,'filename')
    disp(' ');fprintf(1,['Processing dataset from ' indata.filename(1:length(indata.filename)-4) '\n']);
end

%%
UseModel = Fmodel{1};
%Setup Model
Y = UseModel{3}{:};
if ~isempty(UseModel{2})
    %Model has one predictor
    if max(size(UseModel{2}))==1
        X = UseModel{2}{:};
    else
        if istable(UseModel{2})
            X = UseModel{2}{:,:};
        else
            X=[];
            for MPcount = 1 : max(size(UseModel{2}))
                XN = UseModel{2}{MPcount};
                if size(XN,2)>size(XN,1)
                    XN = XN'; %transpose if dimensions do not match
                end
                X = [X XN];
            end
        end
    end
else
    X= [];
end
if size(X,2)>size(X,1)
    X = X'; %transpose if dimensions do not match
end
if size(Y,2)>size(Y,1)
    Y = Y'; %transpose if dimensions do not match
end
if ~isstr(UseModel{3}{:})
    disp(['Using EEG activity as predictor rather than observation...'])
    add_info.EEG_is_predictor=1;
    PName{size(X,2)+1} = 'EEG';
    if FretEEG~=0
        FretEEG=0;
        display('EEG cannot be binned and returned if EEG itself is a regressor... Deactivating.')
    end
end


%Check Designmatrix for NaNs and remove them
L1 = size(X,1);
if ~isempty(find(isnan(X)==1)) || (add_info.EEG_is_predictor && ~isempty(find(isnan(Y)==1)))
    [~,ri]=nanrem(X);
    ri2 = [];
    if add_info.EEG_is_predictor
        [~,ri2]=nanrem(Y);
        Y(unique([ri; ri2]),:)=[];
    end
    Trials(unique([ri; ri2]))=[]; 
    X(unique([ri; ri2]),:)=[];
    if Fdisp > 0
        disp(['Removing ' num2str(L1-size(X,1)) ' rows with NaNs from design matrix.'])
    end
    Warnings(length(Warnings)+1) = {[num2str(L1-size(X,1)) ' NaNs have been removed from Designmatrix.']};
end

%Add information per regressor
N_Pred = size(X,2)+add_info.EEG_is_predictor;
if Level==1
    if Fdisp>0; disp(' ');disp(['Designmatrix includes ' num2str(size(X,1)) ' trials and ' num2str(N_Pred) ' predictors:']); end
else
    if Fdisp>0; disp(' ');disp(['Designmatrix includes ' num2str(size(X,1)) ' subjects and ' num2str(N_Pred) ' predictors:']); end
end

    %%
for Rcount = 1 : size(X,2)
    %Create Dummy Name for Variable if not provided
    if isempty(PredNames)
        PName{Rcount} = {['Predictor ' num2str(Rcount)]};
    else
        PName{Rcount} = PredNames{Rcount};
    end
    %check number of unique entries for regressor
    U = unique(X(:,Rcount));
    if length(U)==2
        PredValues(Rcount) = {U};
        nU = 2; 
        Nentries(Rcount) = {[length(find([X(:,Rcount)]==U(1))) length(find([X(:,Rcount)]==U(2)))]};
        Split = ceil(2* tiedrank(X(:,Rcount)) / sum(~isnan(X(:,Rcount)))); %replace values with 1 and 2
        if Fdisp>0 && ~isempty(PredLables)
            max_stringlength=max(cellfun(@(x) length(x),PredNames));
            spacer_string = ['\t' PName{Rcount}];
            for tabc = 1 : max_stringlength - length(PName{Rcount})
                spacer_string = [spacer_string ' '];
            end
            fprintf([spacer_string ' - Containing ' num2str(Nentries{Rcount}(1)) ' x ' PredLables{Rcount,1} ' (= ' num2str(U(1)) ') and ' num2str(Nentries{Rcount}(2)) ' x ' PredLables{Rcount,2} ' (= ' num2str(U(2)) ')\n'])
        elseif Fdisp>0 && isempty(PredLables)
            fprintf(['Predictor ' num2str(Rcount) ' - Containing ' num2str(Nentries{Rcount}(1)) ' x ' num2str(U(1)) ' and ' num2str(Nentries{Rcount}(2)) ' x ' num2str(U(2)) '\n' ])
        end
    else 
        if FretEEG>0
            nU = FretEEG;
        else
            nU = 3;
        end
        Split = ceil(nU* tiedrank(X(:,Rcount)) / sum(~isnan(X(:,Rcount))));
        for nec = 1 : nU
            NE(nec) = length(find([Split]==nec));
            MeanE(nec) = mean(X(Split==nec,Rcount));
        end
        Nentries(Rcount) = {NE}; clear NE;
        PredValues(Rcount) = {MeanE};
        if Fdisp>0 && ~isempty(PredLables)
            max_stringlength=max(cellfun(@(x) length(x),PredNames));
            spacer_string = ['\t' PName{Rcount}];
            for tabc = 1 : max_stringlength - length(PName{Rcount})
                spacer_string = [spacer_string ' '];
            end
            fprintf([spacer_string ' - Numeric & range ' num2str(min(X(:,Rcount))) ' to ' num2str(max(X(:,Rcount))) ', mean = ' num2str(mean(X(:,Rcount))) '\n'])
        elseif Fdisp>0 && isempty(PredLables)
            fprintf(['Predictor ' num2str(Rcount) '\t - Numeric & range ' num2str(min(X(:,Rcount))) ' to ' num2str(max(X(:,Rcount))) ', mean = ' num2str(mean(X(:,Rcount))) '\n'])
        end
    end
end
if isempty(PredLables)
    for Rcount = 1 : size(X,2)
        PredLables{Rcount,1}  = 'Lable1';
        PredLables{Rcount,2}  = 'Lable2';
    end
end

if Fdisp>0; disp(' ');end
if isempty(IncluReg)
    IncluReg = 1 : length(Nentries);
end
add_info.TotalTrials        = size(X,1);
add_info.RegNames           = PName;
add_info.RegLables          = PredLables;
add_info.RegValues          = PredValues;
add_info.RegNumbers         = Nentries;
add_info.NaNsRemoved        = L1-size(X,1);
add_info.IncluReg           = IncluReg;
% keyboard

%%
for Curr_E = 1 : length(electrodes) %Big loop through all electrodes
    if Curr_E > 1
        TF_Fdisp = 1;
    end 
    if Fdisp == 0
        TF_Fdisp = 0;
    end
    if ~isempty(TF)
        if isfield(TF, 'SecondLevel') && TF.SecondLevel == 1
            RegData = squeeze(indata.data(Curr_E,:,Trials,:));
            ODTF.TF = TF;
            TimeMS = indata.times;
        else
            [RegData, ODTF, ~, TFtimesMs] = STA_Perform_TF(indata, TF, TF_Fdisp, Trials, electrodes(Curr_E),[], TimeInMS); 
            RegData = squeeze(RegData); %has dimensions: time, epoch, frequency
            TimeMS = TFtimesMs;
        end
    else
        RegData = squeeze(indata.data(electrodes(Curr_E),:,Trials));
        TimeMS  = indata.times;
    end
    
    for ReEpochData = 1
        if ~isempty(NewEpoch)
            PadBond = [NewEpoch.boundaries(1) - NewEpoch.padding(1) NewEpoch.boundaries(2) + NewEpoch.padding(end)];
            triallim = find(TimeMS == NewEpoch.Time(1) + PadBond(1)) : find(TimeMS == NewEpoch.Time(1) + PadBond(2));
            RegData_temp = nan(length(triallim), size(RegData,2), size(RegData,3));
            NewTrial=[]; Kap = 0;
            for c = 1 : size(RegData,2)
                if NewEpoch.Time(c)+PadBond(2) > TimeMS(end)
                    disp(['Trial ' num2str(c) ' out of boundaries... removing...'])
                else
                    NewTrial = [NewTrial c];
                    Kap = Kap +1;
                    triallim = find(TimeMS == NewEpoch.Time(c) + PadBond(1)) : find(TimeMS == NewEpoch.Time(c) + PadBond(2));
                    RegData_temp(:,Kap,:) = RegData(triallim, c,:);
                end
            end
            RegData     = RegData_temp; clear RegData_temp
            TimeMS      = PadBond(1) : 1000/indata.srate : PadBond(2);
            TimeInDP    = find(TimeMS == NewEpoch.boundaries(1)) : find(TimeMS == NewEpoch.boundaries(end));
            TimeInMS    = TimeMS(TimeInDP);

            if Fdisp > 0 && Curr_E == 1; disp(' '); disp(['Re-epoching the data: New epochs range from ' num2str(TimeInMS(1)) ' to ' num2str(TimeInMS(end)) 'ms.']); end
            if ~isempty(NewEpoch.baseline)
                if Fdisp > 0; disp(['Using new baseline from ' num2str(NewEpoch.baseline(1)) ' to ' num2str(NewEpoch.baseline(end)) 'ms...']); end
                Bsl     = find(TimeMS == NewEpoch.baseline(1)) : find(TimeMS == NewEpoch.baseline(end));
                BslVal  = repmat(mean(RegData(Bsl,:,:)), 1, size(RegData,1));
                RegData = RegData - BslVal;
            else
                if Fdisp > 0 && Curr_E == 1; disp(['Maintaining baseline from previous epoch...']); end
            end
        else
            NewTrial = 1:size(RegData,2);
            X = X(NewTrial,:);
        end
        add_info.Return_Timewindow_ms = TimeInMS(1:DownSample:end);
        add_info.srate_Hz = indata.srate/DownSample;
    end
    
    %Do this only once on the first iteration through electrodes
    if Curr_E == 1
        if NormaliseOn
            disp(['Normalising design matrix'])
            for N = 1 : size(X,2)
                if length(unique(X(:,N)))==2 %Replace all unique values with �1 if only two are present
                    U=unique(X(:,N));
                    X(find(X(:,N)==U(1)),N)=-1;
                    X(find(X(:,N)==U(2)),N)= 1;
                elseif length(unique(X(:,N)))>2 %all other regressors get normalised
                    X(:,N)=normalise(X(:,N));
                end
            end
        end

        %preallocate
        if ismatrix(RegData) %no TF
            Tstat   = nan(length(electrodes), N_Pred, length(TimeInDP));
            Bstat   = Tstat;
            if RetAll >= 1; R2 = nan(length(electrodes), N_Pred+1, length(TimeInDP)); R2Unique= nan(length(electrodes), N_Pred+1, length(TimeInDP));      end
            if RetAll >= 3; Pstat   = nan(length(electrodes), N_Pred, length(TimeInDP)); Intercept   = nan(length(electrodes), length(TimeInDP));   end
            NormalizationOutput = nan(length(electrodes), 1, 2); %3D = Elect Predictor singleton 2 = mean & SD
        else %preallocate with more dimensions for all frequencies
            Tstat   = nan(length(electrodes), N_Pred, length(TimeInDP), size(RegData,3));
            Bstat   = Tstat;
            if RetAll >= 1; R2      = nan(length(electrodes), N_Pred+1, length(TimeInDP), size(RegData,3)); R2Unique= nan(length(electrodes), N_Pred+1, length(TimeInDP), size(RegData,3)); end   
            if RetAll >= 3; Pstat   = nan(length(electrodes), N_Pred, length(TimeInDP), size(RegData,3)); Intercept = nan(length(electrodes), length(TimeInDP), size(RegData,3));  end 
            NormalizationOutput = nan(length(electrodes), size(RegData,3), 2); %3D = Elect freq 2 = mean & SD
        end
    end
    
    if Fdisp; disp(OutLable{Curr_E}); fprintf(1,'Percent of regression analysis done :    '); end
    for TimeCount = 1 : length(TimeInDP) 
        if Fdisp
            do=floor((TimeCount/(length(TimeInDP)))*100);
            if do < 10
                fprintf(2,'\b%d', do);
            elseif do < 100
                fprintf(2,'\b\b%d', do);
            else
                fprintf(2,'\b\b\b%d', do);
                fprintf(2, '\n');
            end
        end 
        
        tbin         = find(TimeMS == TimeInMS(TimeCount) - FbinsI) : find(TimeMS == TimeInMS(TimeCount) + FbinsI);
        if isempty(tbin)
            display(' ');warning('Time bin after smoothing not found: Input should be milliseconds and needs to align with the sampling rate!'); return;
        end
        if isempty(TF) %time domain
            Nfreq = 1;
            if ~add_info.EEG_is_predictor %EEG itself is observation
                if FbinsI~=0
                    Y = squeeze(mean(RegData( tbin, :)))';
                else
                    Y = RegData( tbin, :)';
                end
            %Add EEG as predictor if it is not set as Y
            elseif add_info.EEG_is_predictor && ~NormaliseOn 
                X(:,N_Pred) = squeeze(mean(RegData( tbin, :)))';
            elseif add_info.EEG_is_predictor && NormaliseOn
                X(:,N_Pred) = normalise(squeeze(mean(RegData( tbin, :))))';
            end
        else
            Nfreq = size(RegData,3);
        end
        
        for Cfreq = 1 : Nfreq
            if ~isempty(TF)
                if ~add_info.EEG_is_predictor
                    Y = squeeze(mean(RegData( tbin, :, Cfreq),1))'; %average over time bin
                end
                %Add EEG as predictor if it is not set as Y
                if add_info.EEG_is_predictor && ~NormaliseOn
                    X(:,N_Pred) = squeeze(mean(RegData(tbin, :, Cfreq),2))';
                elseif add_info.EEG_is_predictor && NormaliseOn
                    X(:,N_Pred) = normalise(squeeze(mean(RegData(tbin, :, Cfreq),1)))';
                end
            end
            %Run that regression!    
            if strcmp(RegMode, 'ols')
                m=fitlm(X, Y(NewTrial));
                cope = m.Coefficients.Estimate;
                tst = m.Coefficients.tStat;
                pval = m.Coefficients.pValue;
            elseif strcmp(RegMode, 'robust')
                warning off
                [cope, stats_rob]   = robustfit(X, Y(NewTrial));
                warning on
                tst=stats_rob.t;
                pval =stats_rob.p; 
            elseif strcmp(RegMode, 'pinv')
                [cope,~,tst] = ols([Y(NewTrial) Y(NewTrial)*3], [ones(length(NewTrial),1) X ]);
                pval = tcdf(-abs(tst), length(NewTrial)-N_Pred-1)*2;
            else
                error(['Unknown argument for RegMode: ' RegMode])
            end
            if RetAll >= 1
                %Simple R squared
                phat = cope(1)+cope(2)*X(:,1);                            %Intercept and first predictor
                if strcmp(RegMode, 'ols')
                    R2(Curr_E,1,TimeCount,Cfreq)= m.Rsquared.Ordinary;
                else
                    R2(Curr_E,2,TimeCount,Cfreq)= corr(Y,phat)^2;
                end
                R2Unique(Curr_E,2,TimeCount,Cfreq)= R2(Curr_E,2,TimeCount); %R2 for only this predictor
                for RR = 1 : length(cope)-2
                    phat = phat + cope(RR+2)*X(:,RR+1);                    %Intercept and other predictors
                    R2(Curr_E,RR+2,TimeCount,Cfreq)= corr(Y,phat)^2;
                    phat2 = cope(1)+cope(RR+2)*X(:,RR+1);                 %Model with this predictor alone
                    R2Unique(Curr_E,RR+2,TimeCount,Cfreq)= corr(Y,phat2)^2;
                end
                R2(Curr_E,1,TimeCount,Cfreq)= corr(Y,phat)^2;
            end
            
            if RetAll >= 3
                Pstat(Curr_E,:,TimeCount,Cfreq)=pval(2:end);
                Intercept(Curr_E,TimeCount,Cfreq)=cope(1);
            end
            Tstat(Curr_E,:,TimeCount,Cfreq)=tst(2:end);
            Bstat(Curr_E,:,TimeCount,Cfreq)=cope(2:end);
        end
    end   
%     keyboard
    %%
    %Return EEG binned by predictors, if set
    if FretEEG 
        X2 = X;
        if isempty(TF) && Curr_E == 1 %pre allocate in first run
            OutEEG = nan(N_Pred, FretEEG, length(electrodes), length(TimeInDP)); %4D predictor, bin, electrode, time
            OutEEG_SD = nan(N_Pred, FretEEG, length(electrodes), length(TimeInDP));
        elseif ~isempty(TF) && Curr_E == 1
            OutEEG = nan(N_Pred, FretEEG, length(electrodes), length(TimeInDP), Nfreq); %5D predictor, bin, electrode, time, frequency
        end
        for Pcount = 1 : N_Pred
            U = X2(:,Pcount);
            U = unique(U(~isnan(U))); %Unique entries for this predictor?
            if length(U)==1
                Warnings(length(Warnings)+1) = {['Regressor: ' PName{Pcount} ' - is constant.']};
                OutValue = [NaN NaN];
            else %this is a parametric regressor
                
                indREM =[];
                if FretEEG == 1 || length(U)==2
                    FC = 2;
                elseif length(U)<=FretEEG
                    FC = length(U);
                else
                    FC = FretEEG;
                end
                
                for BinCount = 1 : FC
                    if length(U)==2
                        Index_Higher = find([X2(:,Pcount)]==U(BinCount));
                        Regr_MeanVal(BinCount) = mean(X2(Index_Higher,Pcount));
                    else %parametric or > 2 levels in this regressors
                        if FretEEG == length(U) %categorial regressor with n entries
                            if BinCount<FC
                                HighCutoff(BinCount) = U(BinCount+1);
                            else
                                HighCutoff(BinCount) = U(BinCount)+1; %this is just because with parametric regressors, the last entry need to be larger than the previous
                            end
                        else
                            HighCutoff(BinCount) = quantile(X2(:,Pcount),(1/FC)*BinCount);
                        end
                        Index_Higher = find(X2(:,Pcount)<HighCutoff(BinCount))';
                        if BinCount > 1
                            Index_Higher = setdiff(Index_Higher, indREM);
                        end
                        indREM = [indREM Index_Higher];
                        Regr_MeanVal(BinCount) = mean(X2(Index_Higher,Pcount));
                        Regr_SeVal(BinCount) = se(X2(Index_Higher,Pcount));
                    end
                    if isempty(TF)
                        OutEEG(Pcount,BinCount,Curr_E,:) = mean(RegData(TimeInDP,Index_Higher),2);
                        OutEEG_SD(Pcount,BinCount,Curr_E,:) = std(RegData(TimeInDP,Index_Higher),0,2);
                    else
                        OutEEG(Pcount,BinCount,Curr_E,:,:) = mean(RegData(TimeInDP,Index_Higher,:),2);
                        if ~isempty(TF) && ~isempty(TF.bands) && length(TF.bands)>2
                            OutEEG_SD(Pcount,BinCount,Curr_E,:,:) = std(RegData(TimeInDP,Index_Higher,:),0,2);
                        else
                            OutEEG_SD(Pcount,BinCount,Curr_E,:,1) = std(RegData(TimeInDP,Index_Higher),0,2);
                        end
                    end
                end
                %%%ERP image%%%
                if ERPimage & isempty(TF) %only in time domain, not TF domain
                    if Fdisp && Curr_E == 1 && Pcount==1
                        ImageSTR =  'Calculating ERP image' ;
                        if ERPimageDim
                            ImageSTR = [ImageSTR ' scaled to ' num2str(ERPimageDim) ' trials'];
                        else
                            ImageSTR = [ImageSTR ' with dimensions equal to regressor length'];
                        end
                        if ERPimageSmooth
                            ImageSTR = [ImageSTR ' and smoothed over ' num2str(ERPimageSmooth) ' vertical trials'];
                        end
                        disp(ImageSTR)
                    end
                    [OutERPimage(Pcount,Curr_E,:,:), TrialValue(:,Pcount)] = erpimage_noplot(RegData(TimeInDP,:), [X2(:,Pcount)]',TimeInMS,[],ERPimageSmooth,ERPimageDim, 'noplot', 'on');
                    if Curr_E==1 && Pcount==1 %pre allocate, unfortunately, only sure about exact size of ERP image after first run :(
                        OutERPimage = nan(N_Pred,length(electrodes),length(TimeInMS),size(OutERPimage,4));
                        OutERPimage(Pcount,Curr_E,:,:) = erpimage_noplot(RegData(TimeInDP,:), [X2(:,Pcount)]',TimeInMS,[],ERPimageSmooth,ERPimageDim, 'noplot', 'on');
                    end
%                     erpimage(RegData(TimeInDP,:), [X2(:,Pcount)]',TimeInMS,[],ERPimageSmooth,ERPimageDim)
                end
                if FC ~= 2
                    OutValue = Regr_MeanVal;
                else
                    OutValue = U;
                end
            end
            OValue{Pcount}=OutValue;
        end
        if Curr_E == length(electrodes)%do at end
            OutEEG(OutEEG==0)=NaN; %replace missing entries with nans instead of zeros
            Odata.('EEG_Values')        = OValue;
            if isempty(TF)
                OutEEG_SD(OutEEG_SD==0)=NaN;
                Odata.EEG_per_regressor = OutEEG(IncluReg,:,:,1:DownSample:end);
                Odata.EEG_SD_per_regressor = OutEEG_SD(IncluReg,:,:,1:DownSample:end);
                if ERPimage
                    if ImageCompression
                        if Fdisp && Curr_E==1; disp('Rounding ERP image to 2 decimals and converting to single precision.');end
                        OutERPimage = single(round(OutERPimage*100)/100);
                    end
                    Odata.ERP_Image     = OutERPimage(:,:,1:DownSample:end,:);
                    Odata.ERP_Image_TV  = TrialValue;
                end
            else
                if Nfreq>1
                    Odata.EEG_per_regressor = OutEEG(IncluReg,:,:,1:DownSample:end,:);
                    if ~isempty(TF) && ~isempty(TF.bands)
                        Odata.EEG_SD_per_regressor = OutEEG_SD(IncluReg,:,:,1:DownSample:end,:);
                    end
                else
                    Odata.EEG_per_regressor = OutEEG(IncluReg,:,:,1:DownSample:end);
                    if ~isempty(TF) && ~isempty(TF.bands)
                        Odata.EEG_SD_per_regressor = OutEEG_SD(IncluReg,:,:,1:DownSample:end);
                    end
                end
            end
        end
        if isempty(TF)
            NormalizationOutput(Curr_E,1,1:2) = [mean(mean(RegData)) mean(std(RegData))]; %mean and mean of SD over all trials per electrode
        else
            NormalizationOutput(Curr_E,:,1:2) = [squeeze(mean(mean(RegData),2)) squeeze(mean(std(RegData),2))];
        end
    end
    if Fdisp > 0; fprintf('\n'); end
end %of electrode loop

if strcmp(TF,'0') | isempty(TF)
    add_info.TF='Time Domain';
else
    add_info.TF=ODTF.TF;
end

add_info.AutocorrMatrix = corrcoef(X);
if ~ReturnInverted && add_info.EEG_is_predictor
    fprintf('\nRemoving all other regressors apart from EEG.')
    RetIndex = size(Tstat,2);
    add_info.RegNames   = add_info.RegNames(end);
    add_info.RegLables  = add_info.RegLables(end,:);
else
    RetIndex = IncluReg;
end

Odata.t_values      = Tstat(:,RetIndex,1:DownSample:end,1:Cfreq);
Odata.b_values      = Bstat(:,RetIndex,1:DownSample:end,1:Cfreq);
Odata.Normalization = NormalizationOutput;
if RetAll >= 1
    Odata.PseudoR2 = R2(:,[1 RetIndex+1],1:DownSample:end,1:Cfreq);
end
if RetAll >= 2
    Odata.UniqueR2 = R2Unique(:,[1 RetIndex+1],1:DownSample:end,1:Cfreq);
end
if RetAll >= 3
    Odata.Robust_p = Pstat(:,RetIndex,1:DownSample:end,1:Cfreq);
    Odata.Intercept = Intercept(:,1:DownSample:end,1:Cfreq);
end

fprintf(2,'\n')
add_info.Warnings = Warnings;
return;































