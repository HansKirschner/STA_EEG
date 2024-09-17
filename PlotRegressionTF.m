function [ hFig, PSavename ] = PlotRegressionTF(D,st,Info,TF)
%%
% keyboard
%%
f_xsize         = 2000;
f_ysize         = 1400;
pHeight         = 0.105;
pWidth          = 0.05;
ShrinkFactor    = 0.8;      %shrinks the plot a little bit
OffSet          = 0.1;      %offest from lower left corner  
NumberOfYTicks  = 4; 
XTickSteps      = 200; %distance between X Tick Labels
%determine small string with short model summary

ti3 = sprintf([strrep(st.ModelName,'_',' ')]);
AbsMinP = 1;
if isempty(st.PlotElec3)
    UElectPlot = st.El2Plot;
else
    UElectPlot = AGF_structfind(Info,'Output_Labels',st.PlotElec3);
end

if length(UElectPlot) > 25
    pHeight         = pHeight/1.5;
    pWidth          = pWidth/1.5;
    NumberOfYTicks  = 2; 
    XTickSteps      = 400; %distance between X Tick Labels
end
XTlable         = st.PlotTime(1):XTickSteps:st.PlotTime(2);
XTposition      = closeval(st.OverallTime,XTlable);

for T = 1 : length(Info.RegNames)
    if ismember(T,st.PlottedRegressors) %this regressor will be plotted
        ti3=[ti3 sprintf(strrep(['\nReg: ' Info.RegNames{T} ' - ' Info.RegLables{T,1} '(' num2str(Info.RegValues{T}(1)) ') & ' Info.RegLables{T,2} '(' num2str(Info.RegValues{T}(2)) ')'  ], '_', ' '))];
    else
        ti3=[ti3 sprintf(strrep(['\nReg (not plotted): ' Info.RegNames{T} ' - ' Info.RegLables{T,1} '(' num2str(Info.RegValues{T}(1)) ') & ' Info.RegLables{T,2} '(' num2str(Info.RegValues{T}(2)) ')'  ], '_', ' '))];
    end
end

if ishandle(1)
    clf; %clear figure 
end
hFig = figure(1);
set(hFig, 'Position', [100 100 f_xsize f_ysize])%
handaxes2 = axes('Position', [OffSet OffSet 1 1]);
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');

text(-pWidth, 0.88,ti3,'HorizontalAlignment','left','VerticalAlignment', 'top', 'FontSize', 11)
clear AllPlots;
AllPlots.Handle=[]; HC=1;
    
if ~st.R2 %plot normal effects
    %Determine maplimits overall
    if st.GroupStats
        [~,~,~,stats]=ttest(D(st.GroupVP,:,st.DataRegNumber,:,:));
        Limi = [-max(max(max(abs(squeeze(stats.tstat))))) max(max(max(abs(squeeze(stats.tstat)))))];
        ti3 = [ti3 sprintf(['\n\nT limits in p = ' num2str(1-tcdf(Limi(2), stats.df(1)))])];
    elseif st.Ninclu~=1
        Limi = [-max(max(max(abs(squeeze(nanmean(D(st.GroupVP,:,st.DataRegNumber,:,:))))))) max(max(max(abs(squeeze(nanmean(D(st.GroupVP,:,st.DataRegNumber,:,:)))))))];
    else
        Limi = [onemax(squeeze(D(st.GroupVP,:,st.DataRegNumber,:,:)),2) onemax(squeeze(D(st.GroupVP,:,st.DataRegNumber,:,:)),1)];
    end
    if st.Level == 1
        ti=['TF Regressor: ' Info.RegNames{st.PlottedRegressors(st.DataRegNumber)} ' - values: ' Info.RegLables{st.PlottedRegressors(st.DataRegNumber),1} '(' num2str(Info.RegValues{st.PlottedRegressors(st.DataRegNumber)}(1)) ') and '...
            Info.RegLables{st.PlottedRegressors(st.DataRegNumber),2} '(' num2str(Info.RegValues{st.PlottedRegressors(st.DataRegNumber)}(2)) ') ' st.UseGroupName]; 
        PSavename = [st.outfold  'Reg Plot ' Info.RegNames{st.PlottedRegressors(st.DataRegNumber)} st.UseGroupName st.AddString];
    else
        ti=['2nd Level TF on ' st.ModelName ' for: ' Info.RegNames{st.PlottedRegressors(st.DataRegNumber)} ' - values: ' Info.RegLables{st.PlottedRegressors(st.DataRegNumber),1} '(' num2str(Info.RegValues{st.PlottedRegressors(st.DataRegNumber)}(1)) ') and '...
            Info.RegLables{st.PlottedRegressors(st.DataRegNumber),2} '(' num2str(Info.RegValues{st.PlottedRegressors(st.DataRegNumber)}(2)) ') ' st.UseGroupName]; 
        PSavename = [st.outfold 'Reg Plot '  Info.RegNames{st.PlottedRegressors(st.DataRegNumber)} ' on ' st.ModelName ' ' st.UseGroupName st.AddString];
        %Overwrite limits
        Limi = [-max(max(max(abs(D(st.GroupVP,:,st.DataRegNumber,:,:))))) max(max(max(abs(D(st.GroupVP,:,st.DataRegNumber,:,:)))))];
    end
    ti = strrep(ti, '_', ' ');
    ti2 = ['Extremes: ' num2str(round2(Limi(1),0.001)) ' to ' num2str(round2(Limi(2),0.001)) ' (' st.WhatIsIt ')'];
    
    %use symmetric map limits (otherwise masking does not work!
    Limi = [-onemax(Limi,1) onemax(Limi,1)];

    if st.MaskPval
        if strcmp(st.MaskPval,'fdr')
            ti2 = [ti2 ' & masked with ' st.MaskPval]; 
        else
            ti2 = [ti2 ' & masked with p < ' num2str(st.MaskPval)]; 
        end
    end
    text(0.5-pWidth, 0.9,ti,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)
    text(0.5-pWidth, 0.87,ti2,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 15)
    
    for E = 1 : length(UElectPlot)
        if isfield(st, 'UseICA') && st.UseICA
            AllX = [1:-(2/length(UElectPlot)):-1];
            ax = AllX(E);
            ay = 0;
        else
            ax =  [st.ChanLoc(UElectPlot(E)).X]; ax(abs(ax)<0.01)=0;
            ay =  [st.ChanLoc(UElectPlot(E)).Y]*ShrinkFactor*-1; ay(abs(ay)<0.01)=0;
        end

        if abs(ax)>5
            ax=ax/100; ay=ay/100; %az=az/100;
        end
        ax = (ax / 2 + 0.5 - pHeight / 2) * ShrinkFactor+OffSet;
        ay = (ay / 2 + 0.5 - pWidth  / 2) * ShrinkFactor+OffSet;
        
        handaxes2 = axes('Position', [ay ax pHeight pWidth]);
        %Determine p-value across subjects
        if st.Level == 1 && st.Ninclu > 1
            [~,PP,~,stats]=ttest(squeeze(D(st.GroupVP,UElectPlot(E),st.DataRegNumber,:,TF.stepnumber:-1:1)));
            PP=squeeze(PP)';
        else
            PP=squeeze(st.p_values(UElectPlot(E),st.RegNumber,:,TF.stepnumber:-1:1))';
        end
        UP = st.MaskPval; p_String = '';
        if strcmpi(st.MaskPval,'fdr')
            UP = PP(:);
            m=length(UP);
            thresh=[1:m]'*st.fdr_q/m;
            sp=sort(UP);
            masked=sp<=thresh;
            IL = find(masked,1,'last');
            if isempty(IL)
                UP=thresh(1); %no significant point
            else
                UP=sp(IL); %set threshold
            end
            p_String = [', crit p = ' num2str(format_pval(UP))];
        end
        
        if st.GroupStats %Plot second level t-stat
            ThePlot = squeeze(stats.tstat)';
        else %Plot average of first level stats
            if st.Level == 1 && st.Ninclu > 1
                ThePlot = [squeeze(nanmean(D(st.GroupVP,UElectPlot(E),st.DataRegNumber,:,TF.stepnumber:-1:1)))]';
            else
                ThePlot = [squeeze(D(1,UElectPlot(E),st.DataRegNumber,:,TF.stepnumber:-1:1))]';
            end
        end

        if st.MaskPval
            ThePlot(PP>UP)=0;
        end
        imagesc(ThePlot, Limi); 
        colormap(st.Fcmap); 
        set(handaxes2, 'Box', 'off')
        handaxes2.XTickLabel = XTlable;
        handaxes2.XTick = XTposition;
        handaxes2.YTick = round(quantile([1:diff(handaxes2.YLim)],NumberOfYTicks));
        handaxes2.YTickLabel = fliplr(round(TF.freq(round(quantile([1:diff(handaxes2.YLim)],NumberOfYTicks)))*10)/10);
        title(st.ChanLoc(UElectPlot(E)).labels);
        if st.PlotPval 
            text(0.01, 0.1,[' min p = ' num2str(format_pval(onemax(PP,2))) p_String],'HorizontalAlignment','left','VerticalAlignment', 'top', 'FontSize', 9);
        end
         
        pause(0.001)
        
        AllPlots(length(AllPlots(1).Handle)+1).Type    = 'TFspectral';
        AllPlots(length(AllPlots(1).Handle)+1).Title   = ['TF spectrum ' Info.RegNames{st.PlottedRegressors(st.RegNumber)}];
        AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
    end

    %create dummy axis with colorbar
    handaxes2 = axes('Position', [0.98 0.90 0.01 pWidth ]);
    imagesc([0 0], Limi);
    set(gca,'Visible','off'); %white image without axes
    colorbar; handaxes2.XTickLabel = []; handaxes2.YTickLabel = []; 
    
    %plot correlation between regressors
    for CorrPlot = 1 : 2 - (st.Ninclu==1 && st.Level == 1)
        handaxes2 = axes('Position', [0.85 0.88 0.09 0.095]);
        if CorrPlot == 1
            AC_plot = squeeze(mean(Info.ACM));
            if st.Ninclu==1 && st.Level == 1
                AC_plot = squeeze(Info.ACM);
            end
            imagesc(AC_plot, [-1 1]); 
            if st.Ninclu==1 && st.Level == 1
                title('r for regressors');
            else
                title('mean r for regressors');
            end
            axis tight; set(gca, 'XTick', [1 : length(Info.RegNames)],  'YTick', [1 : length(Info.RegNames)], 'XTickLabel', AGF_prune_string(Info.RegNames, 13), 'YTickLabel', AGF_prune_string(Info.RegNames, 13))
        elseif CorrPlot > 1 
            [AC_plot] = squeeze(max(Info.ACM));
            [mini] = squeeze(min(Info.ACM));
            AC_plot(abs(AC_plot)<abs(mini)) = mini(abs(AC_plot)<abs(mini));
            imagesc(AC_plot, [-1 1]); 
            title('max / min r for regressors');
            axis tight; set(gca, 'XTick', [1 : length(Info.RegNames)],  'YTick', [], 'XTickLabel', AGF_prune_string(Info.RegNames, 13), 'YTickLabel',{})
        end
        AGF_rotateXLabels(gca,30);
        for col = 1 : length(Info.RegNames)
            for row = 1 : length(Info.RegNames)
                if col~=row
                    text(col,row,num2str(AGF_round2decimals(AC_plot(row,col),2)),'HorizontalAlignment','center','FontSize',6)
                end
            end
        end
    end
    set(hFig,'Units','Points');pos = get(hFig,'Position');set(hFig,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)]);
    saveas (hFig,PSavename,st.printstring);    
else %Plot R2 if present
    if ~isempty(st.R2_values)
        %Determine maplimits overall
        if st.Level == 1 && st.Ninclu ~= 1
            if st.zScoreR2
                Limi = [-max(max(max(abs(squeeze(zscore(mean(st.R2_values(st.GroupVP,:,1,:,:)))))))) max(max(max(abs(squeeze(zscore(mean(st.R2_values(st.GroupVP,:,1,:,:))))))))];
                R2Str = 'z-scored R2';
            else
                Limi = [0 max(max(max(abs(squeeze(mean(st.R2_values(st.GroupVP,:,1,:,:)))))))];
                R2Str = 'actual R2';
            end
        else
            if st.zScoreR2
                Limi = [-max(max(max(abs(squeeze(zscore(st.R2_values(:,1,:,:))))))) max(max(max(abs(squeeze(zscore(st.R2_values(:,1,:,:)))))))];
                R2Str = 'z-scored R2';
            else
                Limi = [0 max(max(max(abs(squeeze(st.R2_values(:,1,:,:))))))];
                R2Str = 'actual R2';
            end
        end
        if st.Level == 1
            ti=['R2 values for the whole model including ' num2str(length(Info.RegNames)) ' regressors ' st.UseGroupName]; 
            PSavename = [st.outfold  ' R2 ' st.UseGroupName st.AddString];
        else
            ti=['2nd level R2 values for the whole model including ' num2str(length(Info.RegNames)) ' regressors ' st.UseGroupName]; 
            PSavename = [st.outfold  ' R2 on ' st.ModelName st.UseGroupName st.AddString];
        end
        ti = strrep(ti, '_', ' ');
        ti2 = ['Maplimits: ' num2str(round2(Limi(1),0.001)) ' to ' num2str(round2(Limi(2),0.001)) ' ' R2Str];
        text(0.5-pWidth, 0.9,ti,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)
        text(0.5-pWidth, 0.87,ti2,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 15)
        for E = 1 : length(UElectPlot)
            if isfield(st, 'UseICA') && st.UseICA
                AllX = [1:-(2/length(UElectPlot)):-1];
                ax = AllX(E);
                ay = 0;
            else
                ax =  [st.ChanLoc(UElectPlot(E)).X]; ax(abs(ax)<0.01)=0;
                ay =  [st.ChanLoc(UElectPlot(E)).Y]*ShrinkFactor*-1; ay(abs(ay)<0.01)=0;
            end
            if abs(ax)>5
                ax=ax/100; ay=ay/100; %az=az/100;
            end
            ax = (ax / 2 + 0.5 - pHeight / 2) * ShrinkFactor+OffSet;
            ay = (ay / 2 + 0.5 - pWidth  / 2) * ShrinkFactor+OffSet;

            handaxes2 = axes('Position', [ay ax pHeight pWidth]);
            set(handaxes2, 'Box', 'off')
            if st.Level == 1 && st.Ninclu~=1
                ThePlot = [squeeze(mean(st.R2_values(st.GroupVP,UElectPlot(E),1,:,TF.stepnumber:-1:1)))]';
            else
                ThePlot = [squeeze(st.R2_values(UElectPlot(E),1,:,TF.stepnumber:-1:1))]';
            end

            if st.zScoreR2
                imagesc(zscore(ThePlot), Limi);
            else
                imagesc(ThePlot, Limi);
            end
            handaxes2.XTickLabel = XTlable;
            handaxes2.XTick = XTposition;
            handaxes2.YTick = round(quantile([1:diff(handaxes2.YLim)],4));
            handaxes2.YTickLabel = fliplr(round(TF.freq(round(quantile([1:diff(handaxes2.YLim)],NumberOfYTicks)))*10)/10);
            title(st.ChanLoc(UElectPlot(E)).labels);
            colormap(st.Fcmap);  
            pause(0.001)
        end
         %create dummy axis with colorbar
        handaxes2 = axes('Position', [0.95 0.93 0.01 pWidth ]);
        imagesc([0 0], Limi);
        set(gca,'Visible','off'); %white image without axes
        colorbar; handaxes2.XTickLabel = []; handaxes2.YTickLabel = []; 
    end
    set(hFig,'Units','Points');pos = get(hFig,'Position');set(hFig,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)]);
    saveas (hFig,PSavename,st.printstring);    
end  

%%
if st.ScalpPower && ~st.R2 %Topographies for collapsed scalp powers for this regressor
    UChans = st.ChanLoc(AGF_structfind(st.ChanLoc,'labels',st.PlotElec2));
    Spalte = length(st.TopoPlotTime) + 3; 
    Zeile = st.BandNumber;
    clear AllPlots;
    AllPlots.Handle=[];

    if ishandle(1)
        clf; %clear figure 
    end 
    hFig = figure(1); set(gcf,'Visible', 'on'); pc=1; HC=1;
    set(hFig, 'Position', [100 100 1800 (Zeile)*200])%
    set(hFig,'PaperPositionMode','Auto')
    figureSize = get(gcf,'Position');
    ti = strrep(['Topographies for collapsed band power for reg ' Info.RegNames{st.PlottedRegressors(st.RegNumber)} '\n in ' st.ModelName st.Gstring st.AddString],'_',' ');
    DatPlot = squeeze(nanmean(D(st.GroupVP,:,st.DataRegNumber,st.TopoPlotTimeDP,:),1));
    UseLim = [-max(abs([min(min(min(DatPlot))) max(max(max(DatPlot)))])) max(abs([min(min(min(DatPlot))) max(max(max(DatPlot)))]))];
    %Determine topoplot maplimits for this band frequency
    for BP = 1 : st.BandNumber
        hand(HC)=subplot(Zeile,Spalte,[pc pc+2]); hold(gca, 'on');HC=HC+1;pc=pc+2;
        PSTR = sprintf(['Collapsed band-power for ' st.SPNames{BP} '\nRanging from ' num2str(round(100*st.usedfreq{BP}(1))/100) ' to ' num2str(round(100*st.usedfreq{BP}(2))/100) ' Hz\nPlotlimits = ' num2str(round(100*UseLim(1))/100) ' to ' num2str(round(100*UseLim(2))/100) ]);
        text(0,0.5,PSTR,'FontSize',12);
        set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');

        for TopoCount = 1 : length(st.TopoPlotTime)
            hand(HC)=subplot(Zeile,Spalte,TopoCount+pc); hold(gca, 'on');
            BandIndex = st.CollapseIndex(BP,1):st.CollapseIndex(BP,2);
            if length(BandIndex) > 1 %this is not only one band
                topoplot(squeeze(mean(DatPlot(:,TopoCount,BandIndex),3)), UChans, 'maplimits', UseLim);colormap(st.Fcmap);   
            else
                topoplot(squeeze(DatPlot(:,TopoCount,BandIndex)), UChans, 'maplimits', UseLim);colormap(st.Fcmap);   
            end
            title([num2str(st.TopoPlotTime(TopoCount)) ' ms']);
            AllPlots(length(AllPlots(1).Handle)+1).Type    = 'TFTopo';
            AllPlots(length(AllPlots(1).Handle)+1).Title   = ['Band' num2str(BP) ' TF Topo ' Info.RegNames{st.PlottedRegressors(st.RegNumber)} ' ' num2str(st.TopoPlotTime(TopoCount)) ' ms'];
            AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
            pause(0.00001)
        end
        pc = pc + TopoCount + 1;
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.99,sprintf(ti),'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)
    set(hFig,'Units','Points');pos = get(hFig,'Position');set(hFig,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)]);
    saveas (hFig,[st.outfold 'Power Topo ' Info.RegNames{st.PlottedRegressors(st.RegNumber)} st.UseGroupName st.AddString],st.printstring)
    
    %Save individual images (may be better for publications etc.)
    if st.IndiviPlot
        %Set order for Items to plot
        %Test if 'Singles' subfolder exists, otherwise create
        if ~exist([st.outfold  '/singles/'])
            mkdir([st.outfold  '/singles/']);
        end
        for c = 1 : length(st.TopoPlotTime)
            if exist('pfig'); pfig = clf; else; pfig = figure; end
            hax_new = copyobj(hand(AllPlots(1).Handle(c)), pfig);
            set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
            colormap(st.Fcmap);
            saveas (pfig,[st.outfold  '/singles/' AllPlots(c).Title], st.printstring);
        end
        close gcf
    end
end




%%
%Plot raw Spectral EEG in a simple table style layout

if st.PlotEEG && isfield(st,'EEG_values') && ~st.R2 && isfield(st,'PlotEEG_Cat_Values')
    NumberOfYTicks  = 4; 
    max_vert_E = 6; %maximum number of electrodes to plot vertically
    needed_plots = ceil(length(UElectPlot)/max_vert_E);
    f_xsize = 1300;
    f_ysize = 150;
    TempE = UElectPlot;
    EEG_Values = st.PlotEEG_Cat_Values{st.RegNumber};
    %map limits over all data at every electrode that should be plotted
    Limi = [onemax(nanmean(st.EEG_values(st.GroupVP,st.RegNumber,:,UElectPlot,:,:)),2) onemax(nanmean(st.EEG_values(st.GroupVP,st.RegNumber,:,UElectPlot,:,:)),1)];

    if st.Level == 1
        ti=['Raw TF per Regressor: ' Info.RegNames{st.PlottedRegressors(st.RegNumber)} ' ' st.UseGroupName]; 
        PSavename = [st.outfold  'raw TF ' Info.RegNames{st.PlottedRegressors(st.RegNumber)} st.UseGroupName st.AddString];
    else
        ti=['2nd Level TF raw TF on ' st.ModelName ' for: ' Info.RegNames{st.PlottedRegressors(st.RegNumber)} ' - values: ' Info.RegLables{st.PlottedRegressors(st.RegNumber),1} '(' num2str(Info.RegValues{st.PlottedRegressors(st.RegNumber)}(1)) ') and '...
            Info.RegLables{st.PlottedRegressors(st.RegNumber),2} '(' num2str(Info.RegValues{st.PlottedRegressors(st.RegNumber)}(2)) ') ' st.UseGroupName]; 
        PSavename = [st.outfold  Info.RegNames{st.PlottedRegressors(st.RegNumber)} ' on ' st.ModelName ' ' st.UseGroupName st.AddString];
    end
    
    
    for c = 1 : needed_plots 
        hFig = clf;
        if length(TempE)>max_vert_E
            TempE2 = TempE(1:max_vert_E);
            TempE(1:max_vert_E) = []; %remove these elcrtodes from temporary electrode array
        else
            TempE2 = TempE(1:end);
        end;
        set(hFig, 'Position', [2615 661 f_xsize f_ysize*length(TempE2)])%
        Zeile = length(TempE2); Spalte = length(EEG_Values)+1;
        k=0;
        for E =  1 : Zeile %count all electrodes in this plot
            for Col = 1 : Spalte
                k=k+1;
                handaxes2=subplot(Zeile, Spalte, k);
                if Col == Spalte
                    if st.PlotEEG_TF_Diff == 1
                        ThePlot = [squeeze(nanmean(st.EEG_values(st.GroupVP,st.RegNumber,Col-1,TempE2(E),:,TF.stepnumber:-1:1)))]' - [squeeze(nanmean(st.EEG_values(st.GroupVP,st.RegNumber,Col-2,TempE2(E),:,TF.stepnumber:-1:1)))]';
                    elseif st.PlotEEG_TF_Diff == 2
                        [~,PP,~,stats] = ttest([squeeze(st.EEG_values(st.GroupVP,st.RegNumber,Col-1,TempE2(E),:,TF.stepnumber:-1:1))], [squeeze(st.EEG_values(st.GroupVP,st.RegNumber,Col-2,TempE2(E),:,TF.stepnumber:-1:1))]);
                        ThePlot = squeeze(stats.tstat)';
                        if st.PlotEEG_TF_maskP
                            ThePlot(squeeze(PP)'>st.PlotEEG_TF_maskP)=0;
                        end
                    end
                    AddTitle = ' - value = difference' ;
                    imagesc(ThePlot, [-onemax(ThePlot,3) onemax(ThePlot,3)]);colorbar
                else
                    if isempty(st.PlotEEG_TF_base) %do not perform additional baseline correction
                        ThePlot = [squeeze(nanmean(st.EEG_values(st.GroupVP,st.RegNumber,Col,TempE2(E),:,TF.stepnumber:-1:1)))]';
                        imagesc(ThePlot,Limi); colorbar
                    else
                        BaseTime = closeval(st.OverallTime,st.PlotEEG_TF_base);
                        ThePlot = ndimsubt(squeeze(st.EEG_values(st.GroupVP,st.RegNumber,Col,TempE2(E),:,TF.stepnumber:-1:1)),2,BaseTime(1):BaseTime(2),'mean','div');
                        imagesc(squeeze(mean(ThePlot))'); colorbar
                    end;
                    AddTitle = [' - value = ' num2str(EEG_Values(Col))];

                end
                set(handaxes2, 'Box', 'off')
                handaxes2.XTickLabel = XTlable;
                handaxes2.XTick = XTposition;
                handaxes2.YTick = round(quantile([1:diff(handaxes2.YLim)],NumberOfYTicks));
                handaxes2.YTickLabel = fliplr(round(TF.freq(round(quantile([1:diff(handaxes2.YLim)],NumberOfYTicks)))*10)/10);
                colormap(st.Fcmap);  
                title(st.ChanLoc(TempE2(E)).labels);
                if E == 1
                    title([st.ChanLoc(TempE2(E)).labels AddTitle])
                else
                    title(st.ChanLoc(TempE2(E)).labels)
                end
                pause(10^-10);
            end
        end
        
        ti = strrep(ti, '_', ' ');
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,ti,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)
        UPSavename = [PSavename ' ' num2str(c) ' of ' num2str(needed_plots)];
        set(hFig,'Units','Points');pos = get(hFig,'Position');set(hFig,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)]);
        saveas (hFig,UPSavename,st.printstring);    
    end
end













