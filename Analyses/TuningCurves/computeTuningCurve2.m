function [Curve, SE, stats, Raw, Trials, outliers] = computeTuningCurve2(Data, StimIndex, varargin)
% Data is numROIs x numTrials
% StimIndex is numTrials x 1

ROIindex = [1,inf];
TrialIndex = [1,inf];

DetermineOutliers = false;
ControlID = 0; % StimID of control trials, or '[]' if no control trial
StimIDs = [];

% saveOut = false;
% saveFile = '';

directory = cd;

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'TrialIndex'
                TrialIndex = varargin{index+1};
                index = index + 2;
            case {'Outliers','DetermineOutliers','outliers'}
                DetermineOutliers = ~DetermineOutliers;
                index = index + 1;
            case 'ControlID'
                ControlID = varargin{index+1};
                index = index + 2;
            case 'StimIDs'
                StimIDs = varargin{index+1};
                index = index + 2;
%             case {'Save', 'save'}
%                 saveOut = true;
%                 index = index + 1;
%             case {'SaveFile', 'saveFile'}
%                 saveFile = varargin{index+1};
%                 index = index + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Data','var') || isempty(Data)
    [Data, p] = uigetfile({'*.rois;*.mat'}, 'Choose ROI file', directory);
    if isnumeric(Data)
        return
    end
    Data = fullfile(p,Data);
end


%% Load data
if ischar(Data)
    load(Data, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = Data;
    end
    Data = ROIdata;
end

% Compute trial means
if isstruct(Data)
    if ~isfield(Data.rois, 'stimMean')
        Data = computeTrialMean(Data);
    end
    Data = gatherROIdata(Data,'stimMean');
end


%% Determine stimuli info
if isempty(StimIDs)
    StimIDs = unique(StimIndex)';
end
numStims = numel(StimIDs);

% Determine order to cycle through stimuli (necessary for t-test)
if ~isempty(ControlID) && any(StimIDs==ControlID)
    controlindex = find(StimIDs==ControlID);                              % locate control ID
    StimIDs = StimIDs([controlindex,setdiff(1:numStims,controlindex)]);   % reorder so control trials are at front
    ControlID = true;
else
    ControlID = false;
end


%% Determine data to analyze

% Determine trials
if islogical(TrialIndex)
    TrialIndex = find(TrialIndex);
elseif TrialIndex(end) == inf
    TrialIndex = cat(2, TrialIndex(1:end-1), TrialIndex(end-1)+1:size(Data,2));
end
Data = Data(:,TrialIndex);
StimIndex = StimIndex(TrialIndex);

% Determine ROIs
if ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:size(Data,1));
end
Data = Data(ROIindex,:);
[numROIs,numTrials] = size(Data);

% Determine outliers for each stimulus
if DetermineOutliers
    fprintf('Determining outliers...');
    outliers = false(numROIs,numTrials);
    for rindex = 1:numROIs
        outliers(rindex,:) = determineOutliers(Data(rindex,:),'GroupID',StimIndex,'type','medianRule');
    end
    Data(outliers) = nan; % remove outliers from dataset
    fprintf('\tComplete\n');
else
    outliers = [];
end

% Determine trials per stim
if isrow(StimIndex)
    StimIndex = StimIndex';
end
if iscolumn(StimIDs)
    StimIDs = StimIDs';
end


%% Calculate average response for each stimulus
fprintf('Computing tuning curves...');

% Compute tuning curves
Trials = arrayfun(@(x) TrialIndex(StimIndex==x), StimIDs, 'UniformOutput',false); % determine trial IDs used for each stim
Raw = arrayfun(@(x) mat2cell(Data(:,StimIndex==x),ones(numROIs,1),nnz(StimIndex==x)), StimIDs, 'UniformOutput',false); % gather data points
Raw = cat(2,Raw{:});
Curve = cellfun(@nanmean, Raw);                         % compute tuning curve
SE = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))), Raw); % compute standard error

% Compute whether each stim is significantly driven
if ControlID
    stats.driven_p = nan(numROIs,numStims);
    stats.driven_p_corr = nan(numROIs,numStims);
    for u = 1:numROIs
        for s = 2:numStims
            [~,stats.driven_p(u,s)] = ttest2(Raw{u,1},Raw{u,s});       % compute t-test btwn current stim and control stim
        end
        [~,~,stats.driven_p_corr(u,2:end)] = fdr_bh(stats.driven_p(u,2:end)); % correct for multiple comparisons
    end
end

% Compute whether each unit is significantly tuned
stats.tuned_p = nan(numROIs,1);
dict = repelem(1:numStims-ControlID, cellfun(@numel, Trials(ControlID+1:end)))'; % creat grouping variable for anova
for u = 1:numROIs
    if all(~isnan(Curve(u,ControlID+1:end)))
        stats.tuned_p(u) = anovan(cat(2,Raw{u,ControlID+1:end})', dict, 'model','full','display','off'); % compute anova on all non-control positions
    end
end
        
fprintf('\tComplete\n');

% 
% %% Save to file
% if saveOut && ~isempty(saveFile)
%     if ~exist(saveFile, 'file')
%         save(saveFile, 'ROIdata', 'outliers', '-mat', '-v7.3');
%     else
%         save(saveFile, 'ROIdata', 'outliers', '-mat', '-append');
%     end
%     fprintf('\tROIdata saved to: %s\n', saveFile);
% end
