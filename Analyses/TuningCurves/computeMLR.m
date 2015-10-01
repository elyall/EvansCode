function [confusionMatrix, stats] = computeMLR(ROIdata, varargin)


FrameIndex = [];
TrialIndex = [1 inf];
ROIindex = [1 inf];
ControlID = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'FrameIndex', 'Frames', 'frames'}
                FrameIndex = varargin{index+1};
                index = index + 2;
            case {'TrialIndex', 'Trials', 'trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
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

if ~exist('ROIdata', 'var')
    [ROIdata, p] = uigetfile({'*.rois;*.mat'}, 'Select ROI file', directory);
    if ~ROIdata
        return
    else
        ROIdata = fullfile(p, ROIdata);
    end
end


%% Load file
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
end


%% Determine parameters
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIindex);

if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:numel(ROIdata.DataInfo.StimID)];
end
numTrials = numel(TrialIndex);

if isempty(FrameIndex)
    FrameIndex = [ROIdata.DataInfo.numFramesBefore+1, ROIdata.DataInfo.numFramesBefore+mode(ROIdata.DataInfo.numStimFrames(TrialIndex))];
end


%% Determine stimuli

% Set control trials to be last stimulus (mnrfit assumes last category as
% reference)
if isempty(ControlID)
    ControlID = min(ROIdata.DataInfo.StimID(TrialIndex));
end
ROIdata.DataInfo.StimID(ROIdata.DataInfo.StimID==ControlID) = max(ROIdata.DataInfo.StimID) + 1;

% Initialize output
[temp,~,StimIndex] = unique(ROIdata.DataInfo.StimID(TrialIndex));
totalStims = numel(temp);


%% Determine number of each type of trial
numStims = nan(totalStims,1);
for sindex = 1:totalStims
    numStims(sindex) = nnz(StimIndex==sindex);
end
minNumStim = min(numStims);


%% Perform analysis
predictions = zeros(numTrials, numROIs);

tic
parfor_progress(numROIs);
parfor rindex = 1:numROIs
    warning('off', 'stats:mnrfit:IterOrEvalLimit');

    % Pull out current ROI's data
    data = mean(ROIdata.rois(ROIindex(rindex)).dFoF(TrialIndex, FrameIndex(1):FrameIndex(2)),2);
    
    % Cycle through testing each trial after training on other trials
    for testIndex = 1:numTrials
        
        % Determine number of trials to use for each stimlus
        stim = StimIndex(testIndex);                    % ID of current trial
        n = min(minNumStim, numStims(stim) - 1);        % minimum number of trials available for single stimuli
        
        % Determine trials to use for each stimulus
        trainIndices = nan(n*totalStims, 1);
        for index = 1:totalStims
            if index~=stim % stimulus is not stimulus being tested
                trainIndices(n*(index-1)+1:n*index) = randsample(find(StimIndex==index), n); % randomly choose n trials for current stimulus
            else % stimulus is one being tested
                temp = find(StimIndex==index);  % find trials for current stimulus
                temp(temp==testIndex) = [];     % remove trial being tested
                trainIndices(n*(index-1)+1:n*index) = randsample(temp, n); % randomly choose n trials for current stimulus
            end
        end
        
        % Compute logistic regression
        [Weights,~,Statistics] = mnrfit(data(trainIndices), StimIndex(trainIndices));
        
        % Validate & generate confusion matrix
        [~,est] = max(mnrval(Weights, data(testIndex), Statistics));
        predictions(testIndex, rindex) = est;
        
    end
    parfor_progress;
end
parfor_progress(0);
toc


%% Build confusion matrices
confusionMatrix = zeros(totalStims, totalStims, numROIs);
for rindex = 1:numROIs
    for sindex = 1:totalStims
        confusionMatrix(sindex, :, rindex) = hist(predictions(StimIndex==sindex,rindex), 1:totalStims);
    end
end


%% Compute stats

% Compute percent correct & selectivity
PCC = nan(numROIs, 1);
Sel = nan(numROIs, 1);
parfor rindex = 1:numROIs
    
    % Normalize confusion matrix
    normCM = bsxfun(@rdivide, confusionMatrix(:,:,rindex), numStims); % numStims == sum(confusionMatrix(:,:,rindex),2)
    
    % Compute PCC
    PCC(rindex) = trace(normCM)/totalStims;
    
    % Compute Selectivity
    [maxV, maxL] = max(diag(normCM));
    Sel(rindex) = log2(maxV*(totalStims-1)/sum(d(setdiff(1:totalStims, maxL)))); % equal to log2(Sel(rindex)/mean(d(~temp2)));
    
end

% Create output
stats.PCC = PCC;
stats.Sel = Sel;

