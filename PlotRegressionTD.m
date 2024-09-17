function [ hFig, PSavename ] = PlotRegressionTD(D,st,Info)
%%
%Image has dimensions: 3 (Info) | n plots | ? 3 Regressor | ? 3 Pval | ? 5 EEG     | ? 5 ERPimage
Spalte = length(st.TopoPlotTime)        +4       + 3*st.PlotRegress + 3*st.PlotPval + 5*st.PlotEEG + 5*st.PlotImage; 
AllPlots.Handle=[];
Zeile = length(st.PlottedRegressors)+1;
if Zeile == 2 %slightly larger plot if only one line of regressors (avoids overlap)
    HeightPerRow = 300;
else
    HeightPerRow = 200;
end
PlotWidth = 2400;
YDirStr = {'normal' 'reverse'};
disp(['Generating ' num2str(length(st.TopoPlotTime)) ' topo plots...'])

if ishandle(1)
    clf; %clear figure 
end

hFig = figure(1); set(gcf,'Visible', 'on'); 
set(hFig, 'Position', [100 100 PlotWidth Zeile*HeightPerRow])%
set(hFig,'PaperPositionMode','Auto')
ti = strrep(['Overview at ' st.PlotElec ' for ' strrep(st.ModelName,'_',' ') st.Gstring st.AddString],'_',' ');
% keyboard
%%
clf
HC = 1; %Handle Counter to export subplots
for c = 1 : Zeile
    disp(['Plotting row ' num2str(c) ' / ' num2str(Zeile)])
    %%%%%%%%%%%R2 Plot and General Information%%%%%%%%%%%%%
    if c == 1 %First row is general R2 of full model
        if isempty(find(isnan(st.R2_values)==1)) && ~isempty(st.R2_values) && st.PlotR2
            %Plot first field with general information
            hand(HC)=subplot(Zeile,Spalte,[1 2]); hold(gca, 'on');HC=HC+1;
            
            if st.Level == 1 && st.Ninclu > 1
                All_Plot_R2 = squeeze(nanmean(st.R2_values(st.GroupVP,:,1,st.TopoPlotTimeDP),1));
            else
                All_Plot_R2 = squeeze(st.R2_values(:,1,st.TopoPlotTimeDP));
            end
            maxR2 = max(max(All_Plot_R2));
            %Determine electrode for R2 plot
            if strcmp(st.PlotElec, 'best')
                [x,~]=find(All_Plot_R2 == maxR2);
                UseElectrode = x;               
            else
                UseElectrode = find(strcmpi({st.ChanLoc.labels},st.PlotElec));
            end
            if isempty(UseElectrode)
                fprintf(['Electrode for plotting (' st.PlotElec ') not found in available electrode set! \nReturned electrodes / ICs are:\n'])
                for plotc = 1 : length(st.ChanLoc)
                    fprintf(['\t' st.ChanLoc(plotc).labels '\n'])
                end
                error('Plotting not possible!')
            elseif length(UseElectrode)>1
                UseElectrode=UseElectrode(1);
            end
            %Make R2 usable for plot
            if st.Level == 1 && st.Ninclu > 1
                All_Plot_R2 = squeeze(nanmean(st.R2_values(st.GroupVP,:,1,st.TopoPlotTimeDP),1));
                All_Plot_R2_full = squeeze(nanmean(st.R2_values(st.GroupVP,UseElectrode,1,:),1));
            else
                All_Plot_R2 = squeeze(st.R2_values(:,1,st.TopoPlotTimeDP));
                All_Plot_R2_full = squeeze(st.R2_values(UseElectrode,1,:));
            end
            if st.zScoreR2
                All_Plot_R2 = zscore(All_Plot_R2);
            else
                All_Plot_R2 = All_Plot_R2 - min(min(All_Plot_R2)); %Subtract minimum R2 (otherwise plot looks bad)                    
            end
            maxR2 = max(max(All_Plot_R2_full));
            minR2 = min(min(All_Plot_R2_full));
            TM = find(All_Plot_R2_full==maxR2);
            
            %Generate string
            if st.Level == 1
                PSTR=[ 'Overrall R2 \nn subjects:  ' num2str(length(st.GroupVP)) '\nn Regressors  ' num2str(length(st.PlottedRegressors)) '\nMin  ' num2str(round2(minR2,0.001)) '\nMax  ' num2str(round2(maxR2,0.001)) ' at ' num2str(st.OverallTime(TM)) ' ms'];
            else
                PSTR=[ 'Overrall R2 \nn subjects:  ' num2str(Info.TotalTrials) '\nn Regressors  ' num2str(length(st.PlottedRegressors)) '\nMin  ' num2str(round2(minR2,0.001)) '\nMax  ' num2str(round2(maxR2,0.001)) ' at ' num2str(st.OverallTime(TM)) ' ms'];
            end
            str = sprintf(PSTR);  
            text(0.1,0.3,str,'FontSize',12);
            set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');
            
            
            
            %Plot R2 topographies
            for pc = 4 : length(st.TopoPlotTimeDP)+3
                hand(HC)=subplot(Zeile,Spalte,pc); hold(gca, 'on');
                set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);
                if isfield(st,'UseICA') && st.UseICA
                    topoplot(st.ICA_Topo{UseElectrode}, st.Original_Chanlocs);colormap(st.Fcmap);        
                    title([strrep(st.ChanLoc(UseElectrode).labels,'_',' ')]);
                else 
                    topoplot(All_Plot_R2(:,pc-3), st.ChanLoc, 'maplimits', [-maxR2 maxR2]);colormap(st.Fcmap); 
                    title([num2str(st.TopoPlotTime(pc-3)) ' ms']);
                end
                AllPlots(length(AllPlots(1).Handle)+1).Type    = 'R2';
                AllPlots(length(AllPlots(1).Handle)+1).Title   = ['R2 ' num2str(st.TopoPlotTime(pc-3)) ' ms'];
                AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
                pause(0.00001);
            end
            
            %Plot R2 time course
            hand(HC)=subplot(Zeile,Spalte,[pc+2 pc+6]); hold(gca, 'on');
            titel=['R2 '];

            if st.Ninclu ~= 1
                This_t = mean(squeeze(st.R2_values(st.GroupVP,UseElectrode,1,:)),1);
                This_SD = st.ConfFactor*(std(squeeze(st.R2_values(st.GroupVP,UseElectrode,1,:)))./sqrt(length(st.GroupVP)-1));
                shadedErrBar(st.OverallTime, This_t, This_SD, st.colorsSingle, st.transp);
            else
                if length(size(st.R2_values)) == 4
                    This_t = squeeze(st.R2_values(1,UseElectrode,1,:));
                else
                    This_t = squeeze(st.R2_values(UseElectrode,1,:));
                end
                This_SD = zeros(length(This_t),1);
                plot(st.OverallTime, This_t, st.colorsSingle);
            end

            if st.Confidence
                title([titel ' shade = ' num2str(st.Confidence*100) '% CI @' st.ChanLoc(UseElectrode).labels]);
            else
                if st.Ninclu == 1 && st.Level == 1
                    title(strrep([titel ' from model  @' st.ChanLoc(UseElectrode).labels],'_',' '));
                else
                    title([titel ' shade = sem']);
                end
            end
            set(gca, 'FontSize', 9); 
            axe = gca;
            axe.XLim = [min(st.OverallTime) max(st.OverallTime)];
            AllPlots(length(AllPlots(1).Handle)+1).Type    = 'R2';
            AllPlots(length(AllPlots(1).Handle)+1).Title   = ['R2 overall time course'];
            AllPlots(1).Handle = [[AllPlots(1).Handle] HC]; HC = HC+1;
        end
        if st.PlotCorr && ~Info.EEG_is_predictor
            if ~exist('pc')
                pc = Spalte -14;
            end
            %plot correlation between regressors
            for CorrPlot = 1 : st.Level*-1+3
                if CorrPlot == 2 && st.Ninclu == 1 && st.Level == 1; break; end
                if CorrPlot == 1
                    hand(HC)=subplot(Zeile,Spalte,[pc+8 pc+10]); hold(gca, 'on');
                    if st.Ninclu == 1 && st.Level == 1
                        title('r for regressors')
                    else
                        title('mean r for regressors')
                    end
                    
                    if st.Ninclu > 1
                        AC_plot = squeeze(mean(Info.ACM(st.GroupVP,:,:)));
                    else
                        AC_plot = squeeze(Info.ACM);
                    end
                else
                    
                    hand(HC)=subplot(Zeile,Spalte,[pc+12 pc+14]); hold(gca, 'on');
                    title('max / min r for regressors');
                    [AC_plot] = squeeze(max(Info.ACM(st.GroupVP,:,:)));
                    [mini] = squeeze(min(Info.ACM(st.GroupVP,:,:)));
                    AC_plot(triu(AC_plot,1)~=0) = mini(triu(mini,1)~=0);
                    
%                     AC_plot=AC_plot';
%                     AC_plot(abs(AC_plot)<abs(mini)) = mini(abs(AC_plot)<abs(mini));
                end
                
                imagesc(AC_plot, [-1 1]); 
                axis tight; set(gca, 'XTick', [1 : length(Info.RegNames)],  'YTick', [1 : length(Info.RegNames)], 'XTickLabel', AGF_prune_string(Info.RegNames, 13), 'YTickLabel', AGF_prune_string( Info.RegNames, 13))
                AGF_rotateXLabels(gca,30);
                for col = 1 : length(Info.RegNames)
                    for row = 1 : 1 : length(Info.RegNames)
                        if col~=row
                            text(col,row,num2str(AGF_round2decimals(AC_plot(row,col),2)),'HorizontalAlignment','center','FontSize',9)
                        end
                    end
                end
                HC=HC+1;
            end
        end
    
    %%%%%%%%%%%%plot single regressor effects%%%%%%%%%%%%%
    else 
        startPOS = (c-1)*Spalte;
        DatPlot = squeeze(nanmean(D(st.GroupVP,:,c-1,st.TopoPlotTimeDP),1));
        if st.Ninclu ~= 1
            [~,DatPlotP,~,stats]=ttest(squeeze(D(st.GroupVP,:,c-1,st.TopoPlotTimeDP)));
            if st.GroupStats
                DatPlot = squeeze(stats.tstat);
            end
        elseif st.PlotPval
            DatPlotP = squeeze(st.p_values(:,c-1,st.TopoPlotTimeDP));
        end
        
        UseLim = [-max(abs([min(min(DatPlot)) max(max(DatPlot))])) max(abs([min(min(DatPlot)) max(max(DatPlot))]))];
        if st.MinMaplimits
            if UseLim(1)>-st.MinMaplimits
                UseLim(1)=-st.MinMaplimits;
            end
            if UseLim(2)<st.MinMaplimits
                UseLim(2)=st.MinMaplimits;
            end
        end
        UseLim = round(UseLim*100)./100;
        [BestElectrode,BestTime]=find([abs(DatPlot)]==max(max(abs(DatPlot(:,:))))); % maximum within the plotted time-frames
        [BestElectrode_PosAllTime,BestTime_PosAllTime]=find([squeeze(nanmean(D(st.GroupVP,:,c-1,:),1))]==max(max(squeeze(nanmean(D(st.GroupVP,:,c-1,:),1))))); % maximum across the whole time
        [BestElectrode_NegAllTime,BestTime_NegAllTime]=find([squeeze(nanmean(D(st.GroupVP,:,c-1,:),1))]==min(min(squeeze(nanmean(D(st.GroupVP,:,c-1,:),1))))); % maximum across the whole time
        SaveBest(c-1)=BestElectrode(1);
        hand(HC)=subplot(Zeile,Spalte,[startPOS+1 startPOS+2]); hold(gca, 'on');HC=HC+1;
        TextHand = gca;
        IRname =  Info.RegNames{st.PlottedRegressors(c-1)};
        if iscell(IRname); IRname = IRname{:}; end %can be user error: if names are cell of cells, correct for this
        
        PSTR=[];
        if ~strcmp(Info.RegLables{st.PlottedRegressors(c-1),2},'numeric') && ~st.EEG_is_predictor
            try PSTR=[ 'Regressor: ' IRname...
                '\n' num2str(Info.RegValues{st.PlottedRegressors(c-1)}(1)) ' = ' num2str(Info.RegLables{st.PlottedRegressors(c-1),1})...
                '\n' num2str(Info.RegValues{st.PlottedRegressors(c-1)}(end)) ' = ' num2str(Info.RegLables{st.PlottedRegressors(c-1),2}), '\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))];
            end
            if ~st.UseICA; PSTR = [PSTR '\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))]; end
        elseif ~strcmp(IRname,'EEG')
            try PSTR=[ 'Regressor: ' IRname...
                    '\nminsplit: ' num2str(round(Info.RegValues{st.PlottedRegressors(c-1)}(1)*100)/100)...
                    '\nmaxsplit: ' num2str(round(Info.RegValues{st.PlottedRegressors(c-1)}(end)*100)/100), '\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))];
            end
            if ~st.UseICA; PSTR = [PSTR '\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))]; end
        else
            try PSTR=[PSTR 'EEG = predictor\nMaplimits ' num2str(UseLim(1)) ' ' num2str(UseLim(2))]
            end
                
        end


        if strcmp(st.PlotElec, 'best')
            UseElectrode = BestElectrode;               
        else
            UseElectrode = find(strcmpi({st.ChanLoc.labels},st.PlotElec));
            if isempty(UseElectrode)
                fprintf(['Electrode for plotting (' st.PlotElec ') not found in available electrode set! \nReturned electrodes / ICs are:\n'])
                for plotc = 1 : length(st.ChanLoc)
                    fprintf(['\t' st.ChanLoc(plotc).labels '\n'])
                end
                error('Plotting not possible!')
            end
        end
        if length(UseElectrode) > 1; UseElectrode=UseElectrode(1); end
        
        %get the p values
        if st.Ninclu ~= 1
            [~,PP]=ttest(squeeze(D(st.GroupVP,UseElectrode,c-1,:)));
        else
            if isfield(st,'p_values')
                PP = squeeze(st.p_values(UseElectrode,c-1,:));
            end
        end
        UP = st.MaskPval;
        if strcmpi(st.MaskPval,'fdr')
            m=length(PP);
            thresh=[1:m]*st.fdr_q/m;
            sp=sort(PP);
            masked=sp<=thresh;
            IL = find(masked,1,'last');
            if isempty(IL)
                UP=thresh(1); %no significant point
            else
                UP=sp(IL); %set threshold
            end
        end
        PSTR = [PSTR '\nCrit p = ' num2str(format_pval(UP))];
        
                %If set, mask all above threshold p-values in plots in white
        if UP<1 & exist('DatPlotP')
            DatPlot(DatPlotP>UP)=0;
        end

        %%%%%%%%%%%%%PLOT TOPOs%%%%%%%%%%%%
        for pc = 4 : length(st.TopoPlotTime)+3
            if isfield(st,'UseICA') && st.UseICA
            else 
                hand(HC)=subplot(Zeile,Spalte,startPOS+pc); hold(gca, 'on');
                topoplot(DatPlot(:,pc-3), st.ChanLoc, 'maplimits', UseLim);colormap(st.Fcmap);        
                title([num2str(st.TopoPlotTime(pc-3)) ' ms']);
                AllPlots(length(AllPlots(1).Handle)+1).Type    = 'Topo';
                AllPlots(length(AllPlots(1).Handle)+1).Title   = ['Topo ' IRname ' ' num2str(st.TopoPlotTime(pc-3)) ' ms'];
                AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1; pause(0.00001);
            end
        end
        CurrPosPlot = pc+1;
        
        if st.PlotRegress
            %Plot regressor time-course
            hand(HC)=subplot(Zeile,Spalte,[startPOS+CurrPosPlot+1 startPOS+CurrPosPlot+2]); hold(gca, 'on');
            CurrPosPlot = CurrPosPlot + 3;
            if st.Ninclu ~= 1
                if st.GroupStats
                    [~,~,~,stats]=ttest(squeeze(D(st.GroupVP,UseElectrode,c-1,:)));
                    This_t = stats.tstat;
                    This_SD = st.ConfFactor*(stats.sd./sqrt(length(st.GroupVP)-1));
                else
                    This_t = mean(squeeze(D(st.GroupVP,UseElectrode,c-1,:)),1);
                    This_SD = st.ConfFactor*(std(squeeze(D(st.GroupVP,UseElectrode,c-1,:)))./sqrt(length(st.GroupVP)-1));
                end
                shadedErrBar(st.OverallTime, This_t, This_SD, st.colorsSingle, st.transp);
            else
                This_t = squeeze(D(st.Ninclu,UseElectrode,c-1,:));
                This_SD = zeros(length(This_t),1);
                plot(st.OverallTime, This_t, st.colorsSingle);
            end
            [M,T] = max(This_t); [Mi,Ti] = min(This_t);
            PSTR = [PSTR '\nMax ' num2str(round2(M, 0.01)) ' at ' num2str(st.OverallTime(T)) ' ms\n'];
            PSTR = [PSTR 'Min ' num2str(round2(Mi, 0.01)) ' at ' num2str(st.OverallTime(Ti)) ' ms\n'];
            set(gca,'YDir',YDirStr{st.InvertYAxis+1}, 'FontSize', 8); 
            axe = gca; axe.XLim = [min(st.OverallTime) max(st.OverallTime)];
            ylabel({[st.WhatIsIt]});
            title(['Regression weight@' strrep(st.ChanLoc(UseElectrode).labels,'_',' ')])
            AllPlots(length(AllPlots(1).Handle)+1).Type    = 'Regressor';
            AllPlots(length(AllPlots(1).Handle)+1).Title   = ['Regressor ' IRname ' time course'];
            AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
            if exist('PP')
               shade_the_back(PP<=UP, [183 183 183]./255, st.OverallTime);
            end
        end

        if st.PlotPval
            %Plot p-value time-course
            hand(HC)=subplot(Zeile,Spalte,[startPOS+CurrPosPlot+1 startPOS+CurrPosPlot+2]); hold(gca, 'on');
            CurrPosPlot = CurrPosPlot + 3;
            
            plot(st.OverallTime, PP, st.colorsPval);
            if st.PlotRegress
               PSTR = [PSTR 'p max = ' format_pval(PP(T)) '\np min = ' format_pval(PP(Ti)) '\n'];
            end
            axe = gca; axe.XLim = [min(st.OverallTime) max(st.OverallTime)];
            set(axe,'yscale','log');
            title(['P-values @' strrep(st.ChanLoc(UseElectrode).labels,'_',' ')]);
            AllPlots(length(AllPlots(1).Handle)+1).Type    = 'Pval';
            AllPlots(length(AllPlots(1).Handle)+1).Title   = ['Pval ' IRname ' time course'];
            AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
        end
        
        %%%%%PLot information string
        %add some more details to the string
%         
        axes(TextHand)
        PSTR=strrep(PSTR, '_', ' ');
        str = sprintf(PSTR);  
        text(0.1,0.5,str,'FontSize',9);
        set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[],'box','off','visible','off');
        
        %Plot EEG over regressor
        if st.PlotEEG
            hand(HC)=subplot(Zeile,Spalte,[startPOS+CurrPosPlot+1 startPOS+CurrPosPlot+4]); hold(gca, 'on');
            CurrPosPlot = CurrPosPlot + 5;
            titel='';

            
            lejenda={};
            try EP = length(Info.RegNumbers{st.PlottedRegressors(c-1)});
            catch EP = length(st.EEG_values(st.PlottedRegressors(c-1)));
            end
            for BinCount = 1:EP
%                             if ~isnan(st.EEG_values(st.GroupVP,c-1,BinCount,UseElectrode,:,:)
                    if st.Level == 1 && st.Ninclu ~= 1
                        This_EEG(BinCount,:)         = nanmean(squeeze(st.EEG_values(st.GroupVP,c-1,BinCount,UseElectrode,:,:)),1);
                        This_EEG_SD(BinCount,:)      = st.ConfFactor*(nanstd(squeeze(st.EEG_values(st.GroupVP,c-1,BinCount,UseElectrode,:,:)))./sqrt(length(st.GroupVP)-1));
                    elseif st.Level == 1 && st.Ninclu == 1
                        This_EEG(BinCount,:)         = squeeze(st.EEG_values(c-1,BinCount,UseElectrode,:));
                        This_EEG_SD(BinCount,:)      = (st.ConfFactor.*squeeze(st.EEG_SD_values(c-1,BinCount,UseElectrode,:)))./sqrt(Info.RegNumbers{st.PlottedRegressors(c-1)}(BinCount)-1);
                    else
                        This_EEG(BinCount,:)         = squeeze(st.EEG_values(c-1,BinCount,UseElectrode,:));
                        This_EEG_SD(BinCount,:)      = (st.ConfFactor.*squeeze(st.EEG_SD_values(c-1,BinCount,UseElectrode,:)))./sqrt(Info.RegNumbers{st.PlottedRegressors(c-1)}(BinCount)-1);
                    end
                    shadedErrBar(st.OverallTime, This_EEG(BinCount,:), This_EEG_SD(BinCount,:), st.colors{BinCount}, st.transp);
                    if Info.RegValues{st.PlottedRegressors(c-1)}(BinCount) > 10
                        lejenda(BinCount)       = {[' m = ' num2str(round(Info.RegValues{st.PlottedRegressors(c-1)}(BinCount)))]};
                    else
                        lejenda(BinCount)       = {[' m = ' num2str(round(100*Info.RegValues{st.PlottedRegressors(c-1)}(BinCount))/100)]};
                    end
%                                 if size(st.EEG_values,3)>2
%                                     if ~isnan(st.EEG_values(1, c-1, 3, 1, 1))
%                                         
%                                     else
%                                         lejenda(1)       = {num2str(Info.RegLables{st.PlottedRegressors(c-1),1})}; lejenda(2)       = {num2str(Info.RegLables{st.PlottedRegressors(c-1),2})};
%                                     end
%                                 else
%                                     lejenda(1)       = {num2str(Info.RegLables{st.PlottedRegressors(c-1),1})}; lejenda(2)       = {num2str(Info.RegLables{st.PlottedRegressors(c-1),2})};
%                                 end

%                             end
            end
            legend(lejenda,'Location','southeast');
            if st.Confidence
                title([titel ' shade = ' num2str(st.Confidence*100) '% CI @' strrep(st.ChanLoc(UseElectrode).labels,'_',' ')]);
            else
                title([titel ' shade = sem']);
            end
            set(gca,'YDir',YDirStr{st.InvertYAxis+1}, 'FontSize', 9); 
            axe = gca;
            axe.XLim = [min(st.OverallTime) max(st.OverallTime)];
            
            if st.Level == 1
                if isfield(Info,'EEG_Factor_Level')
                    axe.YLabel.String=[st.BinPlotTitleStr ' Factor ' Info.EEG_Factor_Level ' '];
                else
                    axe.YLabel.String=[st.BinPlotTitleStr ' '];
                end
            else
                ylabel(['Factor for ' Info.UseValues2ndLevel]);
            end
            AllPlots(length(AllPlots(1).Handle)+1).Type    = 'EEG';
            AllPlots(length(AllPlots(1).Handle)+1).Title   = ['EEG ' IRname ' time course'];
            AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
        end    
        %Plot ERP image for regressor
        if st.ERP_Im_plot
            hand(HC)=subplot(Zeile,Spalte,[startPOS+CurrPosPlot+1 startPOS+CurrPosPlot+4]); hold(gca, 'on');
            CurrPosPlot = CurrPosPlot + 5;    
            SS = abs(st.OverallTime(1)-st.OverallTime(2));
            AddSteps = floor(length(st.OverallTime)/20);
            ERP_I_Time = [st.OverallTime st.OverallTime(end)+SS:SS:st.OverallTime(end)+SS*AddSteps];
            
            if st.Ninclu > 1
                UTW = mean(st.ERPImage_TV(st.GroupVP,:,:));
                z = size(UTW);
                UTW = reshape(UTW,[z(2:end) 1]);
            else
                UTW = st.ERPImage_TV;
            end
            ERP_I_Data = [squeeze(st.ERPImage_values(st.PlottedRegressors(c-1), UseElectrode,:,:)) repmat(norm_and_scale(UTW(st.PlottedRegressors(c-1),:)',onemax(squeeze(st.ERPImage_values(st.PlottedRegressors(c-1), UseElectrode,:,:)),2),onemax(squeeze(st.ERPImage_values(st.PlottedRegressors(c-1), UseElectrode,:,:)),1)),1,AddSteps)];
            imagesc(ERP_I_Time,[1:size(ERP_I_Data,3)],ERP_I_Data);colormap(st.Fcmap);
            axe = gca; 
            axe.XLim = [min(st.OverallTime) max(ERP_I_Time)];
            axe.YLim = [1 size(ERP_I_Data,1)];
            set(gca,'YDir','normal', 'ytick',[axe.YLim(1)+axe.YLim(2)*0.05 axe.YLim(2)-axe.YLim(2)*0.05],'YTickLabel', {'min' 'max', 'xcolor','w','ycolor','w','box','off'});
            ylabel({'trial value'})
            title('ERP image');
            
            AllPlots(length(AllPlots(1).Handle)+1).Type    = 'ERPimage';
            AllPlots(length(AllPlots(1).Handle)+1).Title   = ['ERPimage ' IRname];
            AllPlots(1).Handle = [[AllPlots(1).Handle] HC];HC=HC+1;
        end
    end
end

%Add top title
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,ti,'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 19)

%NEEDS REWORK: CURRENTLY PLOT IS NOT SAVED!
set(hFig,'Units','Points');pos = get(hFig,'Position');set(hFig,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)]);
saveas (hFig,[st.outfold strrep(ti, '_', ' ')],st.printstring);  

%Save individual images (may be better for publications etc.
if st.IndiviPlot
    %Set order for Items to plot
    IPlots = [];
    for c = 1 : length(st.IndivTypes)
        IPlots = [IPlots find(strcmpi({AllPlots.Type},st.IndivTypes{c}))]
    end

    %Test if 'Singles' subfolder exists, otherwise create
    if ~exist([st.outfold  '/singles/'])
        mkdir([st.outfold  '/singles/']);
    end

    for c = 1 : length(IPlots)
        pfig = figure;
        hax_new = copyobj(hand(AllPlots(1).Handle(IPlots(c))), pfig);
        set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
        colormap(st.Fcmap);
        saveas (pfig,[st.outfold  '/singles/' AllPlots(IPlots(c)).Title], st.printstring);
        close(pfig);
    end
end
end