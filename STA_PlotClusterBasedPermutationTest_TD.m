function [] = STA_PlotClusterBasedPermutationTest_TD(allClustInfo,Info2ndLevel,RegressorNr,timesToShow,Electrode)

%% Function that plots statistic map for STA for a specified regressor
% across time (abscissa) and channel (ordinate). Moreover, this function 
% plots Topo plots requested time points and the time course 
% of regression coefficient. 
%
% Input: 
% Info2ndLevel  --> Output from the STA_Plot_Regression function
% RegressorNr   --> which regressor should be plotted
% timesToShow   --> specify the times for the topo plots
% Electrode     --> for which electrode should the coefficient trajectory be
%                   plotted

% set-up colormaps
Fcmap       = load('AGF_cmap.mat');
cbColors=[0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66;...
    0 114 178; 213 94 0; 204 121 167]./256;

% set dynamic color limits
clim = [-1 1]*max(max(abs([allClustInfo{RegressorNr}.clustDetails.threshMap])))*.7;


figure(1);clf;
subplot(3,1,1);hold on
imagesc(Info2ndLevel.time, [], allClustInfo{RegressorNr}.clustDetails.threshMap, [-5, 5]);
colormap(Fcmap.AGF_cmap); 
colorbar
ylabel('channel')
%set(gca,'YTickLabel',Info2ndLevel.ChannelsNames);
ylim([1 length(Info2ndLevel.ChannelsNames)])
xlabel('Time (ms)')
set(gca, 'box', 'off')
set(gca,'clim',clim)
title(string(Info2ndLevel.RegNames(RegressorNr)))

for subp = 1:length(timesToShow)
    subplot(3,length(timesToShow),length(timesToShow)+subp)
    topoplot(allClustInfo{RegressorNr}.clustDetails.threshMap(:, find(Info2ndLevel.time==timesToShow(subp))), Info2ndLevel.Chanloc);
    set(gca,'clim',clim)
    colormap(Fcmap.AGF_cmap);
    colorbar
    title([num2str(timesToShow(subp)) 'ms'])
end


subplot(3,1,3);hold on
plot([Info2ndLevel.time(1), Info2ndLevel.time(end)], [0, 0], '--k')
Timecourse=squeeze(nanmean(Info2ndLevel.regress_values(:,Electrode,RegressorNr,:), 1));
SEM=squeeze(nanstd(Info2ndLevel.regress_values(:,Electrode,RegressorNr,:),1))./sqrt(size(Info2ndLevel.regress_values,1));
H=shadedErrorBar(Info2ndLevel.time,Timecourse',SEM',{'-o', 'color', cbColors(2,:), 'markerfacecolor',cbColors(2,:), 'markerSize', 4, 'lineWidth', 1});
shade_the_back(abs(allClustInfo{RegressorNr}.clustDetails.threshMap(Electrode,:))>0, [183 183 183]./255, Info2ndLevel.time);
ylabel('Coefficient')
xlabel('Time (ms)')
title(['Coefficient Time course @ ' string(Info2ndLevel.ChannelsNames(Electrode))])
set(gca,'box', 'off')

