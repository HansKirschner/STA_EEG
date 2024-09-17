function STA_PlotClusterBasedPermutationTest_TF(data,thresh,Electrode,RegressorNr,numPerm)
%% Function that plots statistic map for STA for a specified regressor in time frequency domain
% 
% Input: 
% Info2ndLevel  --> Output from the STA_Plot_Regression function
% RegressorNr   --> which regressor should be plotted
% Electrode     --> for which electrode should the coefficient trajectory be
%                   plotted
% numPerm       --> Should be at least 1000 for reasonably stable results.

% set up data:
addInfo         = false; %if true this plots some additional info
data.tf         = permute(squeeze(data.regress_values(:,Electrode,RegressorNr,:,:)),[1, 3, 2]);
alpha           = .05/2;
data.frex_Info  = data.frex;
data.timeInfo   = data.times2save;

rng('shuffle')

% initialize null hypothesis maps
permmaps = zeros(numPerm,size(data.tf,2),size(data.tf,3));

         
for permi=1:numPerm
    
    % create a subject length array containing randomly assigned -1's and
    % 1s:
    permArray=ones(size(data.tf, 1), 1);
    permArray(logical(binornd(1, .5, size(data.tf, 1), 1)))=-1;
    
    % multiply each subject timeseries by the -1 or 1 assigned randomly on
    % this permutation
    sz=size(data.tf);
    permMat=repmat(permArray, [1 sz(2:end)]);
    
    permmaps(permi,:,:,:) = squeeze(mean(data.tf.*permMat));
   
end

% plot permutation maps to check distibution - these are expected to center
% around zero and approach gausian
if addInfo
    figure(1);clf;
    subplot(2,1,1)
    histogram(permmaps(:,randperm(length(data.frex_Info),1),randperm(length(data.timeInfo),1)))
    title('PermMap at random Freq and TimePoint')
    subplot(2,1,2)
    histogram(permmaps(:,randperm(length(data.frex_Info),1),randperm(length(data.timeInfo),1)))
    title('PermMap at random Freq and TimePoint')
end

%% compute z- and p-values based on normalized distance to H0 distributions (per pixel)
% p-value
pval = alpha;

% convert p-value to Z value
% if you don't have the stats toolbox, set zval=1.6449;
zval = abs(norminv(pval));

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

if addInfo
    figure(2);clf;
    subplot(2,1,1)
    imagesc(mean_h0)
    axis xy
    title('H0 Map')
    subplot(2,1,2)
    imagesc(std_h0)
    axis xy
    title('STD H0 Map')
end

% now threshold real data...
% first Z-score
zmap = (squeeze(mean(data.tf))-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;


if addInfo

    figure(3), clf

    subplot(221)
    imagesc(data.timeInfo,data.frex_Info,squeeze(mean(data.tf)));
    set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')
    set(gca,'clim',[min(min(squeeze(mean(data.tf)))) max(max(squeeze(mean(data.tf))))])
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('TF map of raw t-values')

    subplot(222)
    imagesc(data.timeInfo,data.frex_Info,zmap)
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('z-map, thresholded no correction')
    set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')

    subplot(223)
    imagesc(data.timeInfo,data.frex_Info,squeeze(mean(data.tf)));
    set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')
    hold on
    contour(data.timeInfo,data.frex_Info,logical(zmap),1,'linecolor','k');
    set(gca,'clim',[min(min(squeeze(mean(data.tf)))) max(max(squeeze(mean(data.tf))))])
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    title('significance regions no correction')

end

%% cluster correction

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,numPerm);

max_cluster_mass = zeros(1,numPerm);


% loop through permutations
for permi = 1:numPerm
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);

        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
            
        tempclustmasses = [];
        for k = 1:numel(islands.PixelIdxList)
        
            tempclustmasses(k) = length(islands.PixelIdxList{k})*abs(nanmean(threshimg(islands.PixelIdxList{k})));

        end

         % store mass of the biggest cluster
        max_cluster_mass(permi) = max(tempclustmasses);
    end
end


% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh      = prctile(max_cluster_sizes,100-(100*thresh));

cluster_thresh_mass = prctile(max_cluster_mass,100-(100*thresh));

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);

zmap_mass = zmap;

for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i})<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end

    if numel(islands.PixelIdxList{i})*abs(nanmean(zmap_mass(islands.PixelIdxList{i})))<cluster_thresh_mass
        zmap_mass(islands.PixelIdxList{i})=0;
    end
end


%%% now some plotting...

figure(4), clf

subplot(321)
imagesc(data.timeInfo,data.frex_Info,squeeze(mean(data.tf)));
set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')
set(gca,'clim',[min(min(squeeze(mean(data.tf)))) max(max(squeeze(mean(data.tf))))])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('TF map of raw values')

subplot(323)
imagesc(data.timeInfo,data.frex_Info,squeeze(mean(data.tf)));
set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')
hold on
contour(data.timeInfo,data.frex_Info,logical(zmap),1,'linecolor','k');
set(gca,'clim',[min(min(squeeze(mean(data.tf)))) max(max(squeeze(mean(data.tf))))])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('significance regions based on cluster size')

subplot(324)
imagesc(data.timeInfo,data.frex_Info,zmap)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('z-map, thresholded based on size')
set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')

subplot(325)
imagesc(data.timeInfo,data.frex_Info,squeeze(mean(data.tf)));
set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')
hold on
contour(data.timeInfo,data.frex_Info,logical(zmap_mass),1,'linecolor','k');
set(gca,'clim',[min(min(squeeze(mean(data.tf)))) max(max(squeeze(mean(data.tf))))])   
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('significance regions based on mass')

subplot(326)
imagesc(data.timeInfo,data.frex_Info,zmap_mass)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('z-map, thresholded based on mass')
set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')

%% final output plot - by default this is based on mass
% set dynamic color limits
clim = [-1 1]*max(max(abs([squeeze(mean(data.tf))])))*.8;
figure(5);clf;
set(0,'defaultAxesFontSize',20);
figure(5)
imagesc(data.timeInfo,data.frex_Info,squeeze(mean(data.tf)));
set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')
hold on
contour(data.timeInfo,data.frex_Info,logical(zmap_mass),1,'linecolor','k');
set(gca,'ylim',[data.frex_Info(1) data.frex_Info(end)],'ydir','norm')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('significance regions based on mass')
colorbar
set(gca,'clim',clim)

return;

