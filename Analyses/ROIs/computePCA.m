function [eigVec,eigVal] = computePCA(ActivityMatrix, varargin)

verbose = false;
ROIindex = [1 inf];
TrialIndex = [1 inf];
StimID = [];

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
[eigVec,eigVal] = eig(cov(ActivityMatrix));


%% Plot trials in PCA space
if verbose
    
    % Project points onto top 3 PCs
    projPts = ActivityMatrix*eigVec(:,end-2:end);
    
    figure;
    if isempty(StimID)
        plot3(projPts(:,1),projPts(:,2),projPts(:,3),'.','MarkerSize',24);
    else
        hold on;
        StimIDs = unique(StimID);
        numStim = numel(StimIDs);
        %stimMean = nan(numStim,3);
        Colors = jet(numStim);
        for sindex=1:numStim
            ci = StimID==StimIDs(sindex);
            plot3(projPts(ci,1),projPts(ci,2),projPts(ci,3),'.','MarkerSize',24,'MarkerEdgeColor',Colors(sindex,:));
            %stimMean(sindex,:) = mean(projPts(ci,:));
        end
        legend(strcat('Stim',{' '},cellstr(num2str(StimIDs))),'Location','BestOutside');
        %plot3(stimMean(:,1),stimMean(:,2),stimMean(:,3),'k','LineWidth',2);
    end
    xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
    grid on
    capVar = diag(eigVal); % variance captured by each eigenvector
    plotVar = sum(capVar(end-2:end)) / sum(capVar);
    title(sprintf('Top 3 PCs that account for %.2f%% of variance',plotVar*100)); 
end

