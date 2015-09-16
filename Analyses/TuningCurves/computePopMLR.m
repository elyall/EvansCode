function [Weights, confusionMatrix, stats] = computePopMLR(ROIdata, varargin)

numKFolds = 10;
numRepeats = 100;

FrameIndex = [];
TrialIndex = [1 inf];
ROIindex = [1 inf];
ControlID = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numKFolds'
                numKFolds = varargin{index+1};
                index = index + 2;
            case 'numRepeats'
                numRepeats = varargin{index+1};
                index = index + 2;
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

if isempty(FrameIndex)
    FrameIndex = [ROIdata.DataInfo.numFramesBefore+1, ROIdata.DataInfo.numFramesBefore+mode(ROIdata.DataInfo.numStimFrames(TrialIndex))];
end


%% Initialize output

% Set control trials to be last stimulus (mnrfit assumes last category as
% reference)
if isempty(ControlID)
    ControlID = min(ROIdata.DataInfo.StimID(TrialIndex));
end
ROIdata.DataInfo.StimID(ROIdata.DataInfo.StimID==ControlID) = max(ROIdata.DataInfo.StimID) + 1;

% Initialize outputs
[StimIDs,~,StimIndex] = unique(ROIdata.DataInfo.StimID(TrialIndex));
numStims = numel(StimIDs);
Weights = zeros(numROIs+1, numStims-1);
confusionMatrix = zeros(numStims, numStims);


%% Format data
Data = ROIs2AvgMat(ROIdata, 'Frames', FrameIndex, 'ROIs', ROIindex);
% [~,order] = sort(Data(:,1), 'descend');
% Data = Data(order,:);


%% Perform analysis
tic
parfor_progress(numRepeats);
parfor bindex = 1:numRepeats
    warning('off', 'stats:mnrfit:IterOrEvalLimit');
    warning('off', 'MATLAB:nearlySingularMatrix');
    
    %% Determine indices for k-means cross validation
    KFoldIndices = nan(numel(TrialIndex), 2);
    for sindex = 1:numStims
        
        % Determine in what trials current stimulus was shown
        currentIndices = find(StimIndex==sindex);
        numTrials = numel(currentIndices);
        
        % Determine number of trials that should exist in each fold
        numPerFold = diff(round(linspace(1,numTrials+1,numKFolds+1)));
        % numPerFold = numPerFold(randperm(numKFolds)); % randomize order of fold sizes (not necessary I believe)
        if any(numPerFold==0)
            warning('Requesting more folds than observations for a given stimulus -> some folds will not have a certain stimulus type.');
        end
        
        % Distribute the trials to the folds
        for kindex = 1:numKFolds
            temp = randsample(numTrials-sum(numPerFold(1:kindex-1)), numPerFold(kindex)); % select random trial(s) for current fold
            KFoldIndices(currentIndices(temp),:) = repmat([kindex, sindex], numPerFold(kindex), 1);
            currentIndices(temp) = [];
        end
        
    end
    
    %% Compute
    currentWeights = zeros(numROIs+1, numStims-1, numKFolds);
    currentMatrix = nan(numStims, numStims, numKFolds);
    for kindex = 1:numKFolds
        
        % Compute logistic regression
        [currentWeights(:,:,kindex),~,stats] = mnrfit(Data(:,TrialIndex(KFoldIndices(:,1)~=kindex))', KFoldIndices(KFoldIndices(:,1)~=kindex,2));
        
        % Validate & generate confusion matrix
        [pihat,~,~] = mnrval(currentWeights(:,:,kindex), Data(:,TrialIndex(KFoldIndices(:,1)==kindex))', stats);
        for sindex = 1:numStims
            currentMatrix(sindex,:,kindex) = mean(pihat(KFoldIndices(KFoldIndices(:,1)==kindex, 2)==sindex,:),1);
        end
        
    end
    
    confusionMatrix = confusionMatrix + mean(currentMatrix,3)/numRepeats;
    Weights = Weights + mean(currentWeights,3)/numRepeats;
    
    parfor_progress;
end
parfor_progress(0);
toc


%% Compute stats
stats = struct();

% Compute percent correct
stats.PCC = trace(confusionMatrix)/numStims;

% Compute selectivity (Elie & Theunissen 2015)
d = diag(confusionMatrix);
[Sel, temp1] = max(d);
temp2 = false(numStims,1);
temp2(temp1) = true;
stats.Sel = log2(Sel*(numStims-1))/sum(d(~temp2)); % equal to log2(Sel(rindex)/mean(d(~temp2)));

