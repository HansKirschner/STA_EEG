function allClustInfo = STA_ClusterBasedPermutationTest(Info2ndLevel,numPerm,neighbors,thresh)
%% funciton that does significance testing while correcting for multiple 
% comparisons via clusterbased permutation testing.
% 
% Input:     
% Info2ndLevel  --> Output from the STA_Plot_Regression function
% numPerm       --> Should be at least 1000 for reasonably stable results.
% neighbors     --> Distance for connection map
%
% Output:
% allClustInfo  --> output of the clusterbased permutation testing

% create a connection mat (N*N matrix of connections between electrodes)
connectionMat = zeros(length(Info2ndLevel.Chanloc));
if  neighbors>0
    for i = 1:length(Info2ndLevel.Chanloc)
        [~,~,~,~,ix] = agf_matchchans(Info2ndLevel.Chanloc,Info2ndLevel.Chanloc(i),'noplot');
        connectionMat(i,ix(1:neighbors))= 1;
    end
end

chanLabels          = Info2ndLevel.Chanloc;
betaMat_forTests    = permute(Info2ndLevel.regress_values, [1,2,4,3]);
xTimes              = Info2ndLevel.time;

% OK, now lets do some stats :)

%loop over regressors
for j = 1:size(betaMat_forTests,4)
    
    % subj x Chann x time
    EEG_dat= betaMat_forTests(:,:,:,j);
   
    % get clusters, cluster sizes, cluster masses for positive clusters (ie
    % p<.05 in a one tailed positive test):
    posClusterInfo=getEEG_clusterSize(EEG_dat, connectionMat, thresh, 'right');
    
    % get clusters, cluster sizes, cluster masses for negative clusters (ie
    % p<.05 in a one tailed negative test):
    negClusterInfo=getEEG_clusterSize(EEG_dat, connectionMat, thresh, 'left');
    
 
    % OK... now we have our statistics of interest... we just need a null
    % distribution to compare them to. Lets run through a loop that:
    % 1) flips the signs of the data for each subject according to a fair coin
    % toss.
    % 2) find the maximum statistic values that you see in the entire dataset
    % from these "permuted" datasets.
    
    for k=1:numPerm
        % create a subject length array containing randomly assigned -1's and
        % 1s:
        permArray=ones(size(EEG_dat, 1), 1);
        permArray(logical(binornd(1, .5, size(EEG_dat, 1), 1)))=-1;
        
        % multiply each subject timeseries by the -1 or 1 assigned randomly on
        % this trial
        sz=size(EEG_dat);
        permMat=repmat(permArray, [1 sz(2:end)]);
        
        % get cluster statistics for permuted dataset:
        permClusterInfo=getEEG_clusterSize(EEG_dat.*permMat, connectionMat, thresh, 'right');
        
        % store the maximum of the statistics, for use in null distribution:
        maxSize(k)=(max(permClusterInfo.clustSizeMap(:)));
        maxWt(k)=(max(permClusterInfo.clustWtMap(:)));
    end
    
    % how big are clusters from permuted (null) data:
    % figure(j)
    % histogram(maxWt, 100)
    % xlabel('Cluster size')

    % For a two tailed test, find the the minimum cluster statistics necessary
    % to beat 97.5% of the null distribution
    % Based on "mass"
    % I like this statistic better!!!
    massTh = prctile(maxWt,  97.5 );
    
    gPos=unique(posClusterInfo.ID_map(posClusterInfo.clustWtMap	>massTh));
    gNeg=unique(negClusterInfo.ID_map(negClusterInfo.clustWtMap	>massTh));
    
    % record details of positive clusters:
    clusterDetails=struct; clusterDetails.isPos =[]; clusterDetails.map=[]; clusterDetails.peakChannel=[]; clusterDetails.peakTime=[];
    for i =1:length(gPos)
        clusterDetails.isPos(i)=true;
        clusterDetails.map{i}=posClusterInfo.ID_map==gPos(i);
        isMax=abs(posClusterInfo.tMap)==  max(abs(posClusterInfo.tMap(posClusterInfo.ID_map==gPos(i))));
        [I,J]=find(isMax&clusterDetails.map{i});
        clusterDetails.peakChannel{i}=chanLabels(I);
        clusterDetails.peakTime(i)=xTimes(J);
    end
    % And negative clusters:
    for i =1:length(gNeg)
        clusterDetails.isPos(end+1)=false;
        clusterDetails.map{end+1}=negClusterInfo.ID_map==gNeg(i);
        isMax=abs(negClusterInfo.tMap)==  max(abs(negClusterInfo.tMap(negClusterInfo.ID_map==gNeg(i))));
        [I,J]=find(isMax&clusterDetails.map{end});
        clusterDetails.peakChannel{end+1}=chanLabels(I);
        clusterDetails.peakTime(end+1)=xTimes(J);
    end
    
    % create thresholded t-map
    threshMap=zeros(size(negClusterInfo.tMap));
    threshMap(posClusterInfo.clustWtMap	>massTh)=posClusterInfo.tMap(posClusterInfo.clustWtMap	>massTh);
    threshMap(negClusterInfo.clustWtMap	>massTh)=negClusterInfo.tMap(negClusterInfo.clustWtMap	>massTh);
    
    sigPosClusters=  unique(posClusterInfo.ID_map(posClusterInfo.clustWtMap>massTh));
    sigNegClusters=  unique(negClusterInfo.ID_map(negClusterInfo.clustWtMap>massTh));
    
    
    % create an array of structures to store permutation info:
    clustInfo.posClusterInfo=posClusterInfo;
    clustInfo.negClusterInfo=negClusterInfo;
    clustInfo.massTh=massTh;
    clustInfo.threshMap=threshMap;
    clustInfo.clustDetails=clustInfo;
    clustInfo.sigPosClusters=sigPosClusters;
    clustInfo.sigNegClusters=sigNegClusters;
    allClustInfo{j}=clustInfo;
    
end
