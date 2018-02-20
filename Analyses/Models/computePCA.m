function [eigVec,eigVal,capVar,projPts] = computePCA(ActivityMatrix, varargin)
% ActivityMatrix is 
verbose = false;
ROIindex = [1 inf];
TrialIndex = [1 inf];
GroupID = [];
Labels = [];
Colors = [];

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'TrialIndex', 'Trials', 'trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
                index = index + 1;
            case {'ID','IDs','Group','idx','StimID'}
                GroupID = varargin{index+1};
                index = index + 2;
            case 'Labels'
                Labels = varargin{index+1};
                index = index + 2;
            case 'Colors'
                Colors = varargin{index+1};
                index = index + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('ActivityMatrix', 'var')
    [ActivityMatrix, p] = uigetfile({'*.rois;*.mat'}, 'Select ROI file', directory);
    if ~ActivityMatrix
        return
    else
        ActivityMatrix = fullfile(p, ActivityMatrix);
    end
end


%% Load file
if ischar(ActivityMatrix)
    ROIFile = ActivityMatrix;
    load(ROIFile, 'ROIdata', '-mat');
    ActivityMatrix = ROIdata; clear ROIdata;
end

%Format data
if isstruct(ActivityMatrix)
    ActivityMatrix = reshapeROIdata(ActivityMatrix,'type','MatMean');
end


%% Determine parameters
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:size(ActivityMatrix,2)];
end

if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:size(ActivityMatrix,1)];
end

% Keep only requested data
ActivityMatrix = ActivityMatrix(TrialIndex,ROIindex);
if ~isempty(GroupID) && length(GroupID) > size(ActivityMatrix,1)
    GroupID = GroupID(TrialIndex);
end


%% Compute PCA
Mean = mean(ActivityMatrix);
ActivityMatrix = ActivityMatrix - Mean;

% [eigVec,eigVal] = eig(cov(ActivityMatrix));
% eigVec = fliplr(eigVec);
% eigVal = rot90(rot90(eigVal));
% capVar = diag(eigVal); 
% capVar = capVar/sum(capVar)*100; % variance captured by each eigenvector
% projPts = ActivityMatrix*eigVec;

[eigVec,projPts,eigVal,~,capVar] = pca(ActivityMatrix);


%% Plot trials in PCA space
if verbose
    
    figure;
    subplot(1,3,1);
    if isempty(GroupID) % plot as if single stimulus
        if isempty(Colors)
            Colors = [0,0,0];
        end
        scatter3(projPts(:,1),projPts(:,2),projPts(:,3),24,Colors,'.');
        
    else % stimuli info provided 
        hold on;
        
        % Determine Stims
        IDs = unique(GroupID);
        if isempty(Labels)
            Labels = strcat('Stim',{' '},cellstr(num2str(IDs)));
        end
        numStim = numel(IDs);
        if isempty(Colors)
            Colors = jet(numStim);
        end
        
        % Plot trials for each stimuli
%         groupMean = nan(numStim,3);
        h = nan(numStim,1);
        for sindex=1:numStim
            ind = GroupID==IDs(sindex);
            h(sindex) = scatter3(projPts(ind,1),projPts(ind,2),projPts(ind,3),12,Colors(sindex,:),'.');
%             groupMean(sindex,:) = mean(projPts(ind,1:3));
%             h(sindex) = scatter3(groupMean(sindex,1),groupMean(sindex,2),groupMean(sindex,3),24,Colors(sindex,:),'*');
        end
        legend(h,Labels,'Location','BestOutside');
        legend boxoff

    end
    xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
    grid on
    title(sprintf('Top 3 PCs that account for %.2f%% of variance',sum(capVar(1:3)))); 
    
    subplot(1,3,2);
    plot(capVar,'o');
    xlabel('principle component');
    ylabel('% variance captured');
    title('Variance captured by each PC');
    
    subplot(1,3,3);
    plot(eigVec(:,1:3));
    legend('PC1','PC2','PC3');
end

