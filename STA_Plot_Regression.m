%%%ENTER HERE CRITERIA
function Info2ndLevel = STA_Plot_Regression( infold, outfold, settings )
%%Function that plots the output of a EEG regression analysis in both the time and time-frequenycc (TF) domain.
%input:
%   infold      = folder with out put of regression analysis
%   outfold     = folder to save plots
%   settings    = struct with fields:
%       .MaskPval       = (Default: 0.05) Set to numeric value, to mask all topo plots with a higher p-value out. Set to 'bonf1' to correct each regressor for every point in time using Bonferroni correction. 
%                         Set to 'bonf2' to correct each regressor for all points in all electrodes (very conservative). Set to 'bonf3' to correct each regressor for every electrode and the total number of plotted regressors (probably too conservative). 
%                         'fdr' uses standard B&H fdr for all p-values (reduce numer of plotted electrodes and regressors for hypothesis testing).    
%       .st.IndiviPlot  = Plot individual images from topographies. When plotting TD data, select whin topos to plot with field "IndivTypes", for TF this is always the projection of power onto scalp. 
%                         Note: This automatically creates a subfolder 'singles' in the output folder. This whole procedure is done because matlab PDF plots in many cases are distorted.
%                         (TODO: Add individual plot for overall TF spectral plot).
%       .st.IndivTypes  = Can be: {'Topo' 'Regressor' 'Pval' 'EEG'} or any combination of these.
%
%
%specific fields for time data:
%       .plotelect       = Electrode used for plotting regression time course and EEG if present besides topo plots ('best' = highest absolute effect; or any electrode name like 'FCz')   

%specific fields for TF data:
%       .PlotElec2          = Cell array of names to plot fot TF scalp projections (e.g. {'FCz' 'Cz' 'C3' 'C4'}, leave empty to use all available electrodes (default).
%       .PlotElec3          = When plotting scalp topographies of power bands, you may likely want to use all electrodes for scalp topographies, but possibly less for full spectral decomposition. In this case, leave
%       .PlotElec2          = {} and enter electrodes to plot in normal projections 
%   Settings for raw power plots (not regression results)
%       .PlotEEG_TF_Diff    = 1 = use difference between conditions, 2 = use t values between conditions (defualt)
%       .PlotEEG_TF_maskP   = masks the condition comparison for TF plots of raw EEG plot with the same p val as all others = []), else enter a numeric p value (e.g., 0.01).
%       .PlotEEG_TF_base    = if set, uses this time window to normalize raw power to (i.e., data is % signal change from this), e.g., [-900 -400]
%                             Note: For differences you may not want an additional baseline, but it is usually recommended for raw power. That's what this setting is for.
%   Settings for bandwise scalp projections
%       .ScalpPower         = if active, projects the frequencies collapsed over '.SPbands' onto the scalp [0/1]. 
%       .SPBands            = Cell of upper and lower cut offs for projection frequencies. Example: {[1 4]; [5 8]; [9 12]; [13 30]}.
%       .SPNames            = Cell of names for frequencies. Example: {'delta' 'theta' 'alpha' beta'}.
%
%
%
%Tipps for TF Data:
%       If you have a very high temporal / spectral resolution or many subjects, data can be very large and cause memory problems. In that case, set '.PlotReg = 1' and then '.PlotReg = 2'... This will only load one regressor per time.
%   'GroupStats'        If first level input, 1 plots second level t values (t-test of all first level effects vs zero), if 0, plots either average within b values or t values (dpending on field 'UseValues')
%%
%requires:  AGF_shadedErrorBar
%           AGF_cmap.mat

%Correct folder name
if ~strcmp(infold(end),'/')
    infold=[infold '/'];
end
if ~strcmp(outfold(end),'/')
    outfold=[outfold '/'];
end
% keyboard
%%
for HideSettings=1
    %Set colors for multicolor plots
    st.colors       = {'g' 'r' 'b' 'c' 'r' 'k'};
    st.colors       = {'m' 'g' 'b' 'c' 'r' 'k'};
    st.colors       = {'g' 'b' 'r' 'c' 'r' 'k'};
    %Set color for regressor values
    st.colorsSingle = 'b';
    st.colorsPval   = 'b';
    
    %Set defaults that are passed to plotting functions
    st.PlotEEG          = 1;                        % Plot EEG timecourse if it is present in output.
    st.transp           = 1;                        % Transparency for shades [1/0]
    st.printstring      = 'pdf';                    % Default fileformat for individual plots.
    st.MinMaplimits     = [];                       % define a minimum for maplimits (e.g. 1)
    st.PlotElec         = 'best';                   % Electrode to plot EEG bins for regressor
    st.Confidence       = 0;                        % CI for shades, 0 = sem, 1 = SD
    st.ERP_Im_plot      = 0;                        % Does not assume an ERP image to be present (overwritten otherwise)
    st.zScoreR2         = 0;                        % Do not zscore R2 value before plotting.
    st.IndiviPlot       = 0;                        % Default: Do not produce individual plots for all subplots
    st.IndivTypes       = {'Topo' 'EEG'};%'Topo' 'Regressor' 'Pval' 'EEG'
    st.BinPlotTitleStr  = 'ERP';
    st.PlotImage        = 1;
    st.PlotRegress      = 1;
    st.PlotPval         = 1; %Default: Plots a small inset minimum p value in TF plots
    st.R2_values        = [];
    st.ScalpPower       = 0;
    st.SPBands          = []; %default: no bands to average
    st.SPBands          = {[1 4]; [5 8]; [9 12]; [13 25]};
    st.SPNames          = {}; %default: no names for bands
    st.SPNames          = {'delta' 'theta' 'alpha' 'beta'};
    st.PlotEEG_TF_Diff  = 2; %1 = use difference between conditions, 2 = use t values between conditions (defualt)
    st.PlotEEG_TF_maskP = [];   %masks the condition comparison for TF plots of raw EEG plot with the same p val as all others
    st.PlotEEG_TF_base  = [];  %if set, uses this time window to normalize raw power to (i.e., data is % signal change from this)
    st.UseICA           = 0;  %assumes normal EEG not component data
    st.ERPImage_values  = []; %ERP Image pre-allocation
    st.ERPImage_TV      = []; %Single trial values for the ERP image
    st.PlotR2           = 1;  %Plot R2 if present in output
    st.PlotCorr         = 1; % plot correlations between regressors
    st.fdr_q            = 0.05; % false discovery rate
    st.InvertYAxis      = 1;  %invert y axis in ERP plots
    
    PlotSteps       = 50;   % plot every x ms
    st.PlotElec2    = {};   % for TF only: These electrodes are plotted
    st.PlotElec3    = {};   % this is a special case: if you want scalp topographies for bandpowers, write in here all electrodes to use for the normal plot
                                % and leace st.PlotElec2 empty
    ExtrapolS       = 100;
    UseValues       = 'b';  % Use 'b' or 't' values from regression output? g = group level t stat  
    MaskPval        = 0.05; % Default: P value amsking is 0.05
    Groups          = {};
    GroupNames      = {''};
    AddString       = ''; 
    TFstring        = '';
    Groupstring     = '';
    SPstring        = '';
    ReduceString    = '';
    GroupStats      = 0;
    NN              = 0;
    FileName        = '';
    Level           = 1; %Default: First level average plot
    normERP         = 0; %Default: Do not normalise ERPs
    TFplotRaw       = 1; %Default: Create a plot for the raw TF values (if they have been returned)
    exclude         = [];   %These participants will be excluded from the analysis
    Confidence      = 0; 
    
    %Read in regression output
    CMap=dir;LP          = mfilename('fullpath');
    if sum(ismember({CMap.name}, 'AGF_cmap.mat')) == 1
        Fcmap       = load('AGF_cmap.mat');
    else
        try
            Fcmap       = load('functions/AGF_regression/AGF_cmap.mat');
        catch
            Fcmap       = load([LP(1:end-length(mfilename)) '/AGF_cmap.mat']);
        end
    end
    st.Fcmap       = Fcmap.AGF_cmap;
    try
        A=load(['data/ChanLocs.mat']);
    catch %load standard channel locations
        if isfield(settings, 'ChanLocs')
            A=load(settings.ChanLocs); %load user provided channel locations
        else
            A=load([LP(1:end-length(mfilename)) '/ChanLocs.mat']); %backup: load a standard file
        end
    end
    if isfield(settings, 'FileName') && ~isempty(settings.FileName)
        F(1).name = settings.FileName;
    else
        F           = dir([infold '*.mat']);
    end
    st.Ninclu      = size(F,1);
    
    %check if files are found?
    if isempty(F)
        disp(['Error: No files found in directory: ' infold]); return
    end
    
    %load first subjects and determine fields of output
    VP              = load([infold F(1).name]);
    FN1             = fieldnames(VP);
    FN2             = fieldnames(VP.(FN1{1}));
    
    %Test if more than one model is saved within the structure
    if size(FN2,1) > 2
        disp('Warning: Found more than one model in input. Models are: ')
        KPMN = 0;
        for PMN = 1 : 2 : size(FN2,1)
            KPMN = KPMN + 1;
            disp(['Model no ' num2str(KPMN) ': ' FN2{PMN}])
        end
        MulM = input('\n\nWhich model do you want to use: ', 's');
        MulM = str2num(MulM);
        if MulM>1
            MulM = MulM * 2 - 1;
        end
        MulI = MulM + 1;
    else
        MulM = 1; MulI = 2;
    end
    
    Info            = VP.(FN1{1}).(FN2{MulI});
    Data            = VP.(FN1{1}).(FN2{MulM});
    if isfield(Data,'Robust_t') %old data input
        TvalFN = 'Robust_t';
    elseif isfield(Data,'Model_1_Robust_t')
        TvalFN = 'Model_1_Robust_t';
    else
        TvalFN = 't_values';
    end
    if isfield(Info, 'Statistic_times_ms') 
        Info.Return_Timewindow_ms = Info.Statistic_times_ms;
    end
    dims                = size(Data.(TvalFN));
    Electrodes          = Info.Output_Labels;
    Elocs               = Info.Output_Electrodes;
    st.EEG_is_predictor = Info.EEG_is_predictor;
    
    AF=fieldnames(A);
    if isfield(Info,'ChannelLocations')
        chanlocs = Info.ChannelLocations;
    else
        chanlocs=A.(AF{1});
    end
    st.ChanLoc = chanlocs(structfind(chanlocs,'labels',Electrodes)); %this is the reduced channel location file that includes only electrodes with data
    if isfield(Info,'ChannelLocations')
        st.Original_Chanlocs = Info.ChannelLocations;
    end
    PlotReg         = [];
    UPlotReg       = 1 : length(Info.RegNames); %Default: plot all regressors
    
    try Info.IncluReg; %%%%%%%%CAN BE REMOVED; FOR COMPATIBILITY ONLY%%%%%%%%
    catch
        Info.IncluReg = 1 : length(Info.RegNames);
    end
    
    %Get input from settings
    if isfield(settings, 'plotsteps');      PlotSteps       = settings.plotsteps;       end
    if isfield(settings, 'confidence');     Confidence      = settings.confidence;      end
    if isfield(settings, 'ExtrapolS');      ExtrapolS       = settings.ExtrapolS;       end
    if isfield(settings, 'MaskPval');       MaskPval        = settings.MaskPval;        end
    if isfield(settings, 'Groups');         Groups          = settings.Groups;          end
    if isfield(settings, 'GroupNames');     GroupNames      = settings.GroupNames;      end
    if isfield(settings, 'UseValues');      UseValues       = settings.UseValues;       end
    if isfield(settings, 'AddString');      AddString       = settings.AddString;       end
    if isfield(settings, 'GroupStats');     GroupStats      = settings.GroupStats;      end
    if isfield(settings, 'PlotReg');        PlotReg         = settings.PlotReg;         end
    if isfield(settings, 'level');          Level           = settings.level;           end
    if isfield(settings, 'normERP');        normERP         = settings.normERP;         end
    if isfield(settings, 'TFplotRaw');      TFplotRaw       = settings.TFplotRaw;       end
    if isfield(settings, 'exclude');        exclude         = settings.exclude;         end
    %Get input from settings ad translate to subfunction settings
    if isfield(settings, 'PlotPval');           st.PlotPval             = settings.PlotPval;        end
    if isfield(settings, 'PlotR2');             st.PlotR2               = settings.PlotR2;          end
    if isfield(settings, 'PlotCorr');           st.PlotCorr             = settings.PlotCorr;        end
    if isfield(settings, 'ScalpPower');         st.ScalpPower           = settings.ScalpPower;      end
    if isfield(settings, 'SPBands');            st.SPBands              = settings.SPBands;         end
    if isfield(settings, 'SPNames');            st.SPNames              = settings.SPNames;         end
    if isfield(settings, 'printstring');        st.printstring          = settings.printstring;     end
    if isfield(settings, 'ERP_Im_plot');        st.ERP_Im_plot          = settings.ERP_Im_plot;     end
    if isfield(settings, 'IndiviPlot');         st.IndiviPlot           = settings.IndiviPlot;      end
    if isfield(settings, 'IndivTypes');         st.IndivTypes           = settings.IndivTypes;      end
    if isfield(settings, 'plotelect2');         st.PlotElec2            = settings.plotelect2;      end
    if isfield(settings, 'plotelect3');         st.PlotElec3            = settings.plotelect3;      end
    if isfield(settings, 'plotImage');          st.plotImage            = settings.plotImage;       end
    if isfield(settings, 'PlotRegress');        st.PlotRegress          = settings.PlotRegress;     end
    if isfield(settings, 'plotPval');           st.plotPval             = settings.plotPval;        end
    if isfield(settings, 'PlotEEG');            st.PlotEEG              = settings.PlotEEG;         end
    if isfield(settings, 'maplimits');          st.MinMaplimits         = settings.maplimits;       end
    if isfield(settings, 'plotelect');          st.PlotElec             = settings.plotelect;       end
    if isfield(settings, 'PlotEEG_TF_Diff');    st.PlotEEG_TF_Diff      = settings.PlotEEG_TF_Diff; end
    if isfield(settings, 'PlotEEG_TF_maskP');   st.PlotEEG_TF_maskP     = settings.PlotEEG_TF_maskP;end
    if isfield(settings, 'PlotEEG_TF_base');    st.PlotEEG_TF_base      = settings.PlotEEG_TF_base; end
    if isfield(settings, 'fdr_q');              st.fdr_q                = settings.fdr_q;           end
    
    
    
    if st.PlotEEG_TF_maskP==2 %has not been set by user
        st.PlotEEG_TF_maskP = st.MaskPval; %default
    end
    
    if ~isempty(Info.IncluReg) %not all regressors are included
        UPlotReg = intersect(Info.IncluReg,UPlotReg);
        ReduceString = [ReduceString 'In the model ' num2str(length(UPlotReg)) ' out of ' num2str(length(Info.RegNames)) ' regressors have been returned.\n'];
    end
    
    if ~isempty(PlotReg)
        [~,UPlotReg] = intersect(UPlotReg,PlotReg);
        UPlotReg = make_wide(UPlotReg);
        if length(UPlotReg) < length(PlotReg)
            disp(['Regressors to be plotted: ' num2str(PlotReg)])
            disp(['Regressors returned in model: ' num2str(Info.IncluReg)])
            error(['Not all regressors that should be plotted have been found in the model!'])
            return;
        end
    else
        PlotReg = UPlotReg;
    end
    ReduceString = [ReduceString num2str(length(UPlotReg)) ' regressors will be plotted.\nThese are:\n'];
    for c = 1 : length(UPlotReg)
        if iscell(Info.RegNames{UPlotReg(c)})
            ReduceString = [ReduceString '\t Reg no ' num2str(PlotReg(c)) ': ' Info.RegNames{UPlotReg(c)}{:} '\n'];
        else
            ReduceString = [ReduceString '\t Reg no ' num2str(PlotReg(c)) ': ' Info.RegNames{UPlotReg(c)} '\n'];
        end
    end
    ReduceString = [ReduceString '\n\n'];
    dims(2) = length(UPlotReg);
    
    %exclude participants if set
    if ~isempty(exclude)
        ReduceString = [ReduceString 'Excluding ' num2str(length(exclude)) ' participants from the analysis.\n'];
        F(exclude)=[];
        st.Ninclu = length(F);
        if length(F) == 1 && Level == 1 %reload initial subject
            VP              = load([infold F(1).name]);
            Info            = VP.(FN1{1}).(FN2{MulI});
            Data            = VP.(FN1{1}).(FN2{MulM});
        end
    end
    
    if GroupStats && Level == 2
        warning('2nd level plotting assumes already group statistics as input. Setting GroupdStats to 0.'); GroupStats = 0;
    end
        
    if st.Ninclu == 1 && Level == 1
        warning('Only one dataset read in but assuming first-level instead of second level output. Is this correct? Otherwise, set settings.level = 2.'); pause(1)
    end
    
    %Test time windows and adapt data
    if ~isfield(settings, 'plottime') %use all available time points
        st.PlotTime = [Info.Return_Timewindow_ms(1) Info.Return_Timewindow_ms(end)];
    else
        st.PlotTime = settings.plottime; %this is the range to be plotted
    end
    %Find these datapoints in the output and generate variables reflecting actual time
    OverallTimeDP = closeval(Info.Return_Timewindow_ms,st.PlotTime(1)) : closeval(Info.Return_Timewindow_ms,st.PlotTime(2));
    st.OverallTime = Info.Return_Timewindow_ms(OverallTimeDP); %this is the exact time
    if isfield(Info,'ICA_Used') && Info.ICA_Used
        st.TopoPlotTime = 1;
        st.UseICA = 1;
        st.ICA_Topo = Info.ICA_Topo;
        ReduceString = [ReduceString '***Using IC activity instead of raw signal as source***\n'];
    elseif length(PlotSteps)>1 %these are indices for times when topoplots should be made (for time domain data only)
        st.TopoPlotTime=PlotSteps;
    else
        st.TopoPlotTime=[st.PlotTime(1) : PlotSteps : st.PlotTime(2)];
    end
    st.TopoPlotTimeDP = closeval(st.OverallTime,st.TopoPlotTime); %Create index for regression datapoints
    if st.PlotTime(1)<Info.Return_Timewindow_ms(1) || st.PlotTime(end)>Info.Return_Timewindow_ms(end)%Test if time windows are okay?
        disp(['Plot time not found in regression output. Stopping. Plottime: ' num2str(st.PlotTime(1)) ' to ' num2str(st.PlotTime(end)) ' and limits of regression: ' num2str(Info.Return_Timewindow_ms(1)) ' to ' num2str(Info.Return_Timewindow_ms(end)) '.']);return;
    end
    
    %%%%%%%%%%%%
    %Pre-allocate arrays, separately for TF or time domain data
    try TF = Info.TF;
    catch
        try TF = Data.TF;
        catch
            TF = 'Time domain';
        end
    end
    Info.ACM  = nan(st.Ninclu,size(Info.AutocorrMatrix,1),size(Info.AutocorrMatrix,1));
    if isstruct(TF)
        if isempty(st.PlotElec2) %default: use all provided electrodes for plot
            st.PlotElec2 = Info.Output_Labels;
        end
        st.El2Plot = structfind(Info,'Output_Labels',st.PlotElec2);
        if length(Info.Output_Labels)==1 %assume that if only one electrode is returned, this one should also be plotted(?)
            st.El2Plot=1;
        end
        if length(dims)>3
            regress_values  = nan(st.Ninclu, length(st.El2Plot), dims(2), length(OverallTimeDP), dims(4));
            if isfield(Data,'PseudoR2')
                st.R2_values = nan(st.Ninclu, length(st.El2Plot), dims(2)+1, length(OverallTimeDP), dims(4));
            end
        else %this happens if TF decomposition was run on only exactly one frequency band
            regress_values  = nan(st.Ninclu, length(st.El2Plot), dims(2), length(OverallTimeDP), 1);
            if isfield(Data,'PseudoR2')
                st.R2_values = nan(st.Ninclu, length(st.El2Plot), dims(2)+1, length(OverallTimeDP), 1);
            end
        end
        TFstring = ['Plotting TF decomposed regression results over ' num2str(length(st.El2Plot)) ' electrodes.\n'];
        TFstring = [TFstring '      ---> Frequency range from ' num2str(TF.frequencies(1)) ' to ' num2str(TF.frequencies(2)) ' Hz with ' num2str(TF.stepnumber) ' ' TF.space ' steps.\n'];
                
        if st.ScalpPower
            if isempty(st.SPBands)
                disp(['To plot scalp power topographies, fequency bands have to be provided in field settings.SPBands']); return;
            end
            st.BandNumber = length(st.SPBands); %number of frequencies to average over
            
            if length(length(st.SPNames))~=st.BandNumber
                disp('Length of names for frequencies to plot does not match number of frequency ranges. Generating new names.')
                st.SPNames = {}; NN = 1;
            end
            
            for CM = 1 : st.BandNumber % find closest matches to these frequencies in the output
                [~,CloseInd1]       = min(abs(TF.freq-st.SPBands{CM}(1)));
                st.usedfreq{CM}(1)     = TF.freq(CloseInd1);
                [~,CloseInd2]       = min(abs(TF.freq-st.SPBands{CM}(2)));
                st.usedfreq{CM}(2)     = TF.freq(CloseInd2);
                st.CollapseIndex(CM,:) = [CloseInd1 CloseInd2];
                if NN %generate new names
                    st.SPNames(CM) = {['Band ' num2str(round(100*st.usedfreq{CM}(1))/100) ' - ' num2str(round(100*st.usedfreq{CM}(2))/100)]};
                end
            end

            SPstring = [SPstring 'Averaging over ' num2str(st.BandNumber) ' frequency ranges for scalp topography plots.\n'];
        end
    else 
        regress_values  = nan(st.Ninclu, dims(1), dims(2), length(OverallTimeDP));
        if isfield(Data,'PseudoR2')
            st.R2_values = nan(st.Ninclu, dims(1), dims(2)+1, length(OverallTimeDP));
        end
        st.El2Plot = 1 : length(Elocs);
    end
    %%%%%%%%%%%%%
    Nsubjects=st.Ninclu;
    if Level==2
        disp(['Model contains 2nd level data of ' num2str(Info.TotalTrials) ' subjects.'])
        if strcmp(Info.UseValues2ndLevel,'t') || strcmp(Info.UseValues2ndLevel,'b')
            st.BinPlotTitleStr = ['first level ' Info.UseValues2ndLevel];
        else
            st.BinPlotTitleStr = ['first level ERP'];
        end
        disp(['2nd level model was run on ' st.BinPlotTitleStr ' data.'])
        Nsubjects=Info.TotalTrials;
    else
         %%%%%%%%%%%%
        %Test if groups are present
        if max([Groups{:}])>st.Ninclu
            disp('Error: Groups contain more values than input folder datasets.');return;
        else
            if ~isempty(Groups)
                Groupstring = ['Comparing ' num2str(length(Groups)) ' groups:\n'];
                for GC = 1 : length(Groups)
                    if ~isempty(GroupNames)
                        Groupstring = [Groupstring '   ' num2str(GC) ' = ' GroupNames{GC} ' with n = ' num2str(length(Groups{GC})) '\n'];
                    else
                        Groupstring = [Groupstring '   ' num2str(GC) ' = Group ' num2str(GC) ' with n = ' num2str(length(Groups{GC})) '\n'];
                        if GC == length(Groups)
                            Groupstring = [Groupstring '\tNote: To specify names for groups, just set value settings.GroupNames.\n'];
                        end
                    end
                end
                Groupstring = [Groupstring '\n'];
            end
        end
    end
    
    %Adjust shades
    if Confidence>0 && Confidence<1
        CString = [num2str(Confidence*100) '% CI'];
        st.ConfFactor = norminv(1+(0.5-(1-Confidence/2)),0,1);  % if confidence is set to calculate CI, scale factor for SE calculation
    elseif Confidence==0
        st.ConfFactor = 1;
        CString = 'SE';
    elseif Confidence==1
        st.ConfFactor = sqrt(Nsubjects-1);
        CString = 'SD';
    end


    
    if normERP == 1
        if ~isfield(Data,'Normalization')
            normERP = 0;
            warning('Cannot perform ERP normalization as output is missing.')
        end
    end
    
    if isfield(Data, 'EEG_per_regressor') && st.PlotEEG == 1
        dims2            = size(VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor);
        if ~isstruct(TF)
            st.EEG_values    = nan(st.Ninclu, length(UPlotReg), dims2(2), dims2(3),length(OverallTimeDP));
        else
            if length(dims2) > 4
                st.EEG_values    = nan(st.Ninclu, length(UPlotReg), dims2(2), dims2(3),length(OverallTimeDP),dims2(end));
            else
                st.EEG_values    = nan(st.Ninclu, length(UPlotReg), dims2(2), dims2(3),length(OverallTimeDP)); 
            end
        end
        if isfield(VP.(FN1{1}).(FN2{MulM}),'EEG_SD_per_regressor')
            st.EEG_SD_values = st.EEG_values;
        end
    else
        if st.PlotEEG == 1;
            warning('No raw EEG signal found in regresison output, yet plot is active? Mistake? Deactivating EEG plot...')
            st.PlotEEG = 0;
        end
    end
    
    if Level == 1
        if strcmp(UseValues,'t')
                st.WhatIsIt = 'average first level t values';
                if GroupStats
                    st.WhatIsIt = 'group level t-stat of within t';
                end
        elseif strcmp(UseValues,'b')
            st.WhatIsIt = 'average first level b values';
            if GroupStats
                st.WhatIsIt = 'group level t-stat of within b';
            end
        end
    else
        st.WhatIsIt = '2nd level t values';
    end
    
    if ~isnumeric(MaskPval)
        PointsFreq = 1;
        if isstruct(TF)
            PointsFreq = dims(end);
        end
        PointsTime = length(OverallTimeDP);
        PointsElec = length(Info.Output_Electrodes);
        PointsReg = length(UPlotReg);
        
        if strcmpi(MaskPval,'bonf1')
            nMC = PointsTime * PointsFreq;
        elseif strcmpi(MaskPval,'bonf2')
            nMC = PointsTime * PointsFreq * PointsElec;
        elseif strcmpi(MaskPval,'bonf3')
            nMC = PointsTime * PointsFreq * PointsElec * PointsReg;
        elseif strcmpi(MaskPval,'fdr')
            %FDR requires sample p values to be used, we just leave the value
        else
            warning(['Unknown setting for p value masking (' MaskPval ') - setting default of 0.05.'])
            MaskPval = 0.05; nMC=[];
        end
        if exist('nMC') && ~isempty(nMC)
            MaskPval = 0.05 / nMC;
            Groupstring = [Groupstring 'Plots masked for ' num2str(nMC) ' multiple comparisons using Bonferroni correction.\n'];
        end
    elseif MaskPval < 0 || MaskPval > 1
        Groupstring = [Groupstring 'P value masking not possible because value: ' num2str(MaskPval) '.\n'];
        MaskPval = 1;
    end
    
    %Display information on prompt about input and model
    fprintf('\n\n************Starting Plot Regression Script...***********\n')
    if Level == 1
        disp(['Found ' num2str(st.Ninclu) ' files in folder ' infold '.'])
    end
    fprintf(Groupstring)
    disp(['Regresison model ' strrep(FN2{MulM},'_',' ') ' includes ' num2str(length(Info.Output_Electrodes)) ' electrodes and ' num2str(min(Info.Return_Timewindow_ms)) ' until ' num2str(max(Info.Return_Timewindow_ms)) 'ms.'])
    disp(['Number of regressors is ' num2str(length(Info.RegNames)) ' and time window for plotting ' num2str(min(st.PlotTime)) ' to ' num2str(max(st.PlotTime)) ' ms.'])
    fprintf(TFstring)
    fprintf([ReduceString])
    fprintf(SPstring)
    
    
    if isfield(Data, 'ERP_Image') && st.ERP_Im_plot
        dims3       = size(VP.(FN1{1}).(FN2{MulM}).ERP_Image);
        if ExtrapolS < dims3(4) %more data available than extrapolation should be --> do not extrapolate!
            ExtrapolS = dims3(4);
        end
        
        if prod([st.Ninclu length(UPlotReg) dims3(2) ExtrapolS dims3(3)])>10^9 %VERY large array, use slower cumulative averaging instead
            %Dimensions are:     Regressor | Electrode | Image Y | Image X
            st.ERPImage_values = nan(length(UPlotReg), dims3(2), ExtrapolS, dims3(3));
            st.ERPImage_TV =  nan(st.Ninclu, length(UPlotReg), ExtrapolS); %trial wise values are 3D: n VP, n reg, length of image
            disp('********')
            disp(['Warning: ERP images are very large, using cumulative averaging (can be slow)!'])
            disp('********')
            
        else
            if Level==2 || (st.Ninclu == 1 && Level == 1)
                st.ERPImage_values= nan(length(UPlotReg), dims3(2), ExtrapolS, dims3(3));
                st.ERPImage_TV =  nan(length(UPlotReg), ExtrapolS); %trial wise values are 2D: n reg, length of image
            else
                %Dimensions are: VP | Regressor | Electrode | Image Y | Image X
                st.ERPImage_values= nan(st.Ninclu, length(UPlotReg), dims3(2), ExtrapolS, dims3(3));
                st.ERPImage_TV =  nan(st.Ninclu, length(UPlotReg), ExtrapolS); %trial wise values are 3D: n VP, n reg, length of image
            end
        end
        
        disp(['Will extrapolate ERP image(s) to ' num2str(ExtrapolS) ' datapoints.' ])
    elseif st.ERP_Im_plot
        disp('ERP image plot data not found. Deactivating...')
        st.ERP_Im_plot = 0;
    end
end
clear Data %remove dummy Data from first subject
% keyboard
%%
%Load data
if Level==2 || (st.Ninclu == 1 && Level == 1)
    if isfield(VP.(FN1{1}).(FN2{MulM}),'Robust_p') && st.PlotPval
        st.p_values        = VP.(FN1{1}).(FN2{MulM}).Robust_p(st.El2Plot,UPlotReg,OverallTimeDP,:);
    elseif  st.PlotPval
       error('No p-values from analysis returned! Set ''retall'' argument for regression call to 3 to activate this!'); return; 
    end
    if isfield(VP.(FN1{1}).(FN2{MulM}), 'PseudoR2') 
        st.R2_values    = VP.(FN1{1}).(FN2{MulM}).PseudoR2(st.El2Plot,1,OverallTimeDP,:);
    end
    
    st.Ninclu = 1;
    Groups = [];
    
    if ~isstruct(TF)
        regress_values  = nan([1, size(VP.(FN1{1}).(FN2{MulM}).b_values(st.El2Plot,:,OverallTimeDP))]);
        if st.PlotEEG == 1
            st.EEG_values = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,st.El2Plot,OverallTimeDP);  %4D --> Regressor | bins | electrode | time
        end
        if isfield(VP.(FN1{1}).(FN2{MulM}),'EEG_SD_per_regressor') && st.PlotEEG == 1
            st.EEG_SD_values = VP.(FN1{1}).(FN2{MulM}).EEG_SD_per_regressor(UPlotReg,:,st.El2Plot,OverallTimeDP);
        end
    else
        if length(size(VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor))>4 %more than one frequency
            st.EEG_values = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,:,OverallTimeDP,:); %5D --> Regressor | bins | electrode | time | freq
            if isfield(VP.(FN1{1}).(FN2{MulM}), 'EEG_SD_per_regressor')
                st.EEG_SD_values = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,:,OverallTimeDP,:); 
            end
        else
            st.EEG_values = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,:,OverallTimeDP); %4D --> Regressor | bins | electrode | time
            if isfield(VP.(FN1{1}).(FN2{MulM}), 'EEG_SD_per_regressor')
                st.EEG_SD_values = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,:,OverallTimeDP); 
            end
        end
    end
    if Level == 2
        AddString = [AddString ' for ' Info.ModelName]; 
    end
end

if st.Ninclu>1; fprintf('Loading data from: '); end
for p = 1 : st.Ninclu
    %Get fitting INFO file and display it
    VPname = [F(p).name(1:end-4) '  '];
    if p > 1 && p~=st.Ninclu
        fprintf(2,[RemoveString '%s'], VPname);
    elseif p == st.Ninclu && st.Ninclu>1
        fprintf(2,[RemoveString '%s'], VPname);
        fprintf('\n\n')
    elseif st.Ninclu>1
        fprintf(2,VPname);
    end
    RemoveString = [VPname VPname]; RemoveString(1:2:length(RemoveString))='\';RemoveString(2:2:length(RemoveString))='b';
    VP = load([infold F(p).name]);
    %extract correlation information between regressors
    Info.ACM(p,:,:) = VP.(FN1{1}).(FN2{MulM+1}).AutocorrMatrix;
    
    %load actual data into arrays
    if isstruct(TF)
        %TF plot discards electrodes that are not included in the plot (save memory)
        if strcmp(UseValues,'t')
            regress_values(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).(TvalFN)(st.El2Plot,UPlotReg,OverallTimeDP,:);
        elseif strcmp(UseValues,'b')
            regress_values(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).b_values(st.El2Plot,UPlotReg,OverallTimeDP,:);
        else
            disp(['Unrecognized value for parameter UseValues: ' UseValues]);return;
        end
        
        if ~isempty(st.R2_values) && Level == 1 && st.Ninclu ~= 1
            st.R2_values(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).PseudoR2(st.El2Plot,[1 UPlotReg],OverallTimeDP,:);
        end
        
    else %Load Timedomain Data
        %Save to grandtotal matrix
        if strcmp(UseValues,'t')
            regress_values(p,:,:,:) = VP.(FN1{1}).(FN2{MulM}).(TvalFN)(:,UPlotReg,OverallTimeDP);
        elseif strcmp(UseValues,'b')
            regress_values(p,:,:,:) = VP.(FN1{1}).(FN2{MulM}).b_values(:,UPlotReg,OverallTimeDP);
        else
            disp(['Unrecognized value for parameter UseValues: ' UseValues]);return;
        end
        
        if ~isempty(st.R2_values)  && Level == 1 && st.Ninclu ~= 1
            st.R2_values(p,:,:,:) = VP.(FN1{1}).(FN2{MulM}).PseudoR2(:,[1 UPlotReg],OverallTimeDP);
        end
        
        if exist('dims3') && st.ERP_Im_plot
            for Rc = 1 : length(UPlotReg)
                for Ec = 1 : size(Electrodes,2)
                    GridInt = griddedInterpolant([squeeze(VP.(FN1{1}).(FN2{MulM}).ERP_Image(UPlotReg(Rc),Ec,:,:))]');
                    vec     = AGF_exact_veclength(1, size([squeeze(VP.(FN1{1}).(FN2{MulM}).ERP_Image(UPlotReg(Rc),Ec,:,:))]',1), ExtrapolS);
                    [Xq,Yq] = ndgrid(vec,1:dims3(3));
                    if prod(size(st.ERPImage_values),2) > 10^9 %use cumulative averaging
                        if p==1
                            st.ERPImage_values(Rc,Ec,:,:) = GridInt(Xq,Yq);
                        else
                            st.ERPImage_values(Rc,Ec,:,:) = (squeeze(st.ERPImage_values(UPlotReg(Rc),Ec,:,:)).*(p-1) + GridInt(Xq,Yq)) / p;
                        end
                    else
                        if Level==2 || (st.Ninclu == 1 && Level == 1)
                            st.ERPImage_values(Rc,Ec,:,:) = GridInt(Xq,Yq);
                        else
                            st.ERPImage_values(p,Rc,Ec,:,:) = GridInt(Xq,Yq);
                         end
                    end
                    if (Level==2 || (st.Ninclu == 1 && Level == 1)) && Ec==1 %ERP Image values are the same for every electrode
                        try %this should work if all participants have the minimum trial number
                            st.ERPImage_TV(Rc,:) = [VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(1,Rc); VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(:,Rc); VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(end,Rc)];
                        catch %interpolate
                            V = [VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(1,Rc); VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(:,Rc); VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(end,Rc)];
                            st.ERPImage_TV(Rc,:) = V(2:end-1);%trial wise values are 2D: n reg, length of image
                        end
                    elseif Ec == 1 %multi subject data
                        st.ERPImage_TV(p,Rc,:) = [VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(1,Rc); VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(:,Rc); VP.(FN1{1}).(FN2{MulM}).ERP_Image_TV(end,Rc)];
                        
                    end
                end
            end
        end
    end
    %Load raw EEG / TF data
    if st.PlotEEG
        if Level == 1 && st.Ninclu ~= 1
            if ~isstruct(TF) %time domain
                st.EEG_values(p,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,:,OverallTimeDP); %5D --> VP | Regressor | bins | electrode | time
                if isfield(VP.(FN1{1}).(FN2{MulM}),'EEG_SD_per_regressor') 
                    st.EEG_SD_values = VP.(FN1{1}).(FN2{MulM}).EEG_SD_per_regressor(UPlotReg,:,:,OverallTimeDP);
                end
            else
                st.EEG_values(p,:,:,:,:,:) = VP.(FN1{1}).(FN2{MulM}).EEG_per_regressor(UPlotReg,:,:,OverallTimeDP,:); %6D --> VP | Regressor | bins | electrode | time | freq
            end
            CountCat = 0;
            for c = 1 : dims(2)
                Ireg = UPlotReg(c);
                dimsEEG = size([VP.(FN1{1}).(FN2{MulM}).EEG_Values{Ireg}]);
                if dimsEEG(1)<dimsEEG(2) %parametric regressor
                    CountCat = CountCat +1;
                    for c2 = 1 : dimsEEG(2)
                        CatReg(p,CountCat,c2) = VP.(FN1{1}).(FN2{MulM}).EEG_Values{Ireg}(c2); % cannot easily be pre-allocated as the size depends on the maximum number of splits in EEG data
                    end
                end
                if p == length(st.Ninclu) %average values for parametric regressors
                    if dimsEEG(1)<dimsEEG(2)
                        st.PlotEEG_Cat_Values{c} = squeeze(mean(CatReg(:,CountCat,:),1));
                    else
                        st.PlotEEG_Cat_Values{c} = VP.(FN1{1}).(FN2{MulM}).EEG_Values{Ireg};
                    end
                end
            end
        end
    end
end
%%
%Average ERP images if not done before...
if prod(size(st.ERPImage_values),2) < 10^9 && Level == 1 && st.Ninclu ~= 1
    z = size(mean(st.ERPImage_values,1));
    st.ERPImage_values = reshape(mean(st.ERPImage_values,1),[z(2:end) 1]); %cannot squeeze with only one electrode, need to reshape to remove first dimension, but none of the others
end
%Adjust for inverted Regression
if strcmp(Info.RegNames{end},'EEG')
    Info.RegLables{end+1}={'numeric'; 'numeric'};
end
if ~exist(outfold)
    mkdir(outfold);
end
% keyboard
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%START%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check if multiple groups are specified? 
GL = length(Groups);
if GL == 0 
    GL = 1;
else
    disp(['Plotting for ' num2str(GL) ' different groups.'])
end

close all; %close other plots
figure;
for GroupCount = 1 : GL
    %Addgroupname if missing
    st.Gstring = ''; st.UseGroupName='';
    if isempty(Groups)
        GroupVP = 1 : st.Ninclu;
    else
        GroupVP = Groups{GroupCount};
        try 
            if isempty(GroupNames{GroupCount})
                st.UseGroupName = ['Group' num2str(GroupCount)];
                st.Gstring = [' in group ' st.UseGroupName];
            elseif ~strcmp(GroupNames{GroupCount}(1),' ')
                st.UseGroupName = [' ' GroupNames{GroupCount}];
                st.Gstring = [' in group ' st.UseGroupName];
            else
                st.UseGroupName = GroupNames{GroupCount};
                st.Gstring = [' in group ' st.UseGroupName];
            end
        catch
            st.UseGroupName = ['Group' num2str(GroupCount)];
            st.Gstring = '';
        end
    end
    
    clear AllPlots
    %Time Frequency Plots
    st.ModelName            = FN2{MulM};
    st.GroupStats           = GroupStats;
    st.PlottedRegressors    = PlotReg; %the plot function (still) requires regressors in the order they are returned
    st.GroupVP              = GroupVP;
    st.Level                = Level;
    st.outfold              = outfold;
    st.AddString            = AddString;
    st.MaskPval             = MaskPval;
    st.UseValues            = UseValues;
    
    if isstruct(TF) && isempty(TF.bands) %Call time frequency plotter for full spectral data
        st.R2        = 0; %do not plot R2 here
        for R = 1 : length(UPlotReg)
            %st.RegNumber = UPlotReg(R);
            st.RegNumber        = UPlotReg(R); %index to all regressors
            st.DataRegNumber    = R; %index to included regressors in data
            PlotRegressionTF(regress_values,st,Info,TF);
        end
        if isempty(find(isnan(st.R2_values))) && st.PlotR2 %Plot R2
            st.R2               = 1; %activate R2 plot
            PlotRegressionTF(regress_values,st,Info,TF);
        end
    elseif isstruct(TF)
        st.InvertYAxis      = 0;
        for BandC = 1 : length(TF.bands)
            st.AddString = [st.AddString ' ' TF.bands{BandC}];
            if length(length(TF.bands)) == 1
                PlotRegressionTD(regress_values,st,Info);
            else
                PlotRegressionTD(regress_values(:,:,:,:,BandC),st,Info);
            end
        end
    else %ERP-like Plots
        PlotRegressionTD(regress_values,st,Info);
    end
end

%% extract info for other 2nd level analysis
Info2ndLevel=struct;

Info2ndLevel.regress_values   = regress_values;
Info2ndLevel.time             = st.OverallTime;
Info2ndLevel.Chanloc          = st.ChanLoc;
Info2ndLevel.originalChanloc  = st.Original_Chanlocs;
Info2ndLevel.RegNames         = Info.RegNames;
Info2ndLevel.ChannelsNr       = Info.Output_Electrodes;
Info2ndLevel.ChannelsNames    = Info.Output_Labels;

if isstruct(TF) && isempty(TF.bands)

    Info2ndLevel.frex       = TF.freq;
    Info2ndLevel.times2save = Info.Return_Timewindow_ms;

end

return;
%%