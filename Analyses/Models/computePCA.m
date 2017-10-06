function [eigVec,eigVal,capVar] = computePCA(ActivityMatrix, varargin)

verbose = false;
ROIindex = [1 inf];
TrialIndex = [1 inf];
StimID = [];
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
            case 'StimID'
                StimID = varargin{index+1};
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
if ~isempty(StimID) && length(StimID) > size(ActivityMatrix,1)
    StimID = StimID(TrialIndex);
end


%% Compute PCA
Mean = mean(ActivityMatrix);
ActivityMatrix = ActivityMatrix - Mean;
[eigVec,eigVal] = eig(cov(ActivityMatrix));
capVar = diag(eigVal); 
capVar = capVar/sum(capVar); % variance captured by each eigenvector


%% Plot trials in PCA space
if verbose
    
    % Project points onto top 3 PCs
    projPts = ActivityMatrix*eigVec(:,end-2:end);
    
    figure;
    if isempty(StimID) % plot as if single stimulus
        if isempty(Colors)
            Colors = [0,0,0];
        end
        scatter3(projPts(:,1),projPts(:,2),projPts(:,3),24,Colors);
        
    else % stimuli info provided 
        hold on;
        
        % Determine Stims
        StimIDs = unique(StimID);
        if isempty(Labels)
            Labels = strcat('Stim',{' '},cellstr(num2str(StimIDs)));
        end
        numStim = numel(StimIDs);
        if isempty(Colors)
            Colors = jet(numStim);
        end
        
        % Plot trials for each stimuli
        stimMean = nan(numStim,3);
        h = nan(numStim,1);
        for sindex=1:numStim
            ind = StimID==StimIDs(sindex);
            h(sindex) = scatter3(projPts(ind,1),projPts(ind,2),projPts(ind,3),12,Colors(sindex,:),'filled');
            stimMean(sindex,:) = mean(projPts(ind,:));
%             h(sindex) = scatter3(stimMean(sindex,1),stimMean(sindex,2),stimMean(sindex,3),24,Colors(sindex,:),'*');
        end
        legend(h,Labels,'Location','BestOutside');
        legend boxoff

    end
    xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
    grid on
    plotVar = sum(capVar(end-2:end));
    title(sprintf('Top 3 PCs that account for %.2f%% of variance',plotVar*100)); 
end

