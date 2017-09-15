function [Weights, confusionMatrix, stats] = computePopMLR(Data, StimIndex, varargin)
% Data is numROIs by numTrials

numKFolds = 5;
numRepeats = 1;
ControlID = [];

% For loading only
FrameIndex = [];
TrialIndex = [];
ROIindex = [];      % ROIindex is a list of ROI indices corresponding to which ROIs to pull data from
FileIndex = [];     % FileIndex is a list of file indices corresponding to which files the ROIs to pull from are in

directory = cd;

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
            case 'ControlID'
                ControlID = varargin{index+1};
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
            case 'FileIndex'
                FileIndex = varargin{index+1};
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

if ~exist('Data', 'var') || isempty(Data)
    [Data, p] = uigetfile({'*.rois;*.mat'},'Select ROI file',directory,'MultiSelect','on');
    if ~Data
        return
    end
    Data = fullfile(p, Data);
end


%% Load data
if ~isnumeric(Data)
    [Data, StimIndex] = ROIs2AvgMat(Data,...
        'ROIindex',     ROIindex,...
        'FileIndex',    FileIndex,...
        'FrameIndex',   FrameIndex,...
        'TrialIndex',   TrialIndex);
else
    if ~isempty(ROIindex)
        Data = Data(ROIindex,:);
    end
    if ~isempty(TrialIndex)
        Data = Data(:,TrialIndex);
        StimIndex = StimIndex(TrialIndex);
    end
end
numROIs = size(Data,1);


%% Initialize output

% Set control trials to be last stimulus (mnrfit assumes last category as reference)
if isempty(ControlID)
    ControlID = min(StimIndex);
end
StimIndex(StimIndex==ControlID) = max(StimIndex) + 1;
[StimIDs,~,StimIndex] = unique(StimIndex);

% Initialize outputs
numStims = numel(StimIDs);
Weights = zeros(numROIs+1, numStims-1);
confusionMatrix = zeros(numStims, numStims);


%% Perform analysis
parfor_progress(numRepeats*numKFolds);
for n = 1:numRepeats
%     warning('off', 'stats:mnrfit:IterOrEvalLimit');
%     warning('off', 'MATLAB:nearlySingularMatrix');
    
    % Determine indices for k-means cross validation
    KFoldIndices = crossvalind('KFold', StimIndex, numKFolds);
    
    % Compute MLR
    currentWeights = zeros(numROIs+1, numStims-1, numKFolds);
    currentMatrix = nan(numStims, numStims, numKFolds);
    for kindex = 1:numKFolds
        
        % Compute logistic regression
        [currentWeights(:,:,kindex),~,stats] = mnrfit(Data(:,KFoldIndices~=kindex)', StimIndex(KFoldIndices~=kindex));
        
        % Validate & generate confusion matrix
        [pihat,~,~] = mnrval(currentWeights(:,:,kindex), Data(:,KFoldIndices==kindex)', stats);
        for sindex = 1:numStims
            currentMatrix(sindex,:,kindex) = mean(pihat(StimIndex(KFoldIndices==kindex)==sindex,:),1);
        end
        
        parfor_progress;
    end
    
    confusionMatrix = confusionMatrix*(n-1)/n + mean(currentMatrix,3)/n; % compute running average of confusion matrix
    Weights = Weights*(n-1)/n + mean(currentWeights,3)/n;                % compute running average of weights
    
end
parfor_progress(0);


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

