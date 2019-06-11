function [Weights, confusionMatrix, stats] = computePopMLR(Data, StimIndex, varargin)
% Data is numROIs by numTrials

numKFolds = 5;
numRepeats = 1;
% ControlID = [];

% For loading only
FrameIndex = [];
TrialIndex = [];
ROIindex = [];      % ROIindex is a list of ROI indices corresponding to which ROIs to pull data from
FileIndex = [];     % FileIndex is a list of file indices corresponding to which files the ROIs to pull from are in
verbose = false;

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
%             case 'ControlID'
%                 ControlID = varargin{index+1};
%                 index = index + 2;
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
            case 'verbose'
                verbose = varargin{index+1};
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
        StimIndex = StimIndex(TrialIndex,:);
    end
end
[numROIs,numTrials] = size(Data);

% % Set control trials to be last stimulus (mnrfit assumes last category as reference)
% if isempty(ControlID)
%     ControlID = min(StimIndex);
% end
% StimIndex(StimIndex==ControlID) = max(StimIndex) + 1;

% Determine stimuli indices
numS = size(StimIndex,2);
[StimIDs,~,Index] = unique(StimIndex,'rows');
numStims = size(StimIDs,1);


%% Perform analysis

% Initialize output
Weights = zeros(numROIs+1, numStims-1, numRepeats);
confusionMatrix = zeros(numStims, numStims, numRepeats);

if verbose; parfor_progress(numRepeats*numKFolds); end
parfor n = 1:numRepeats
%     warning('off', 'stats:mnrfit:IterOrEvalLimit');
%     warning('off', 'MATLAB:nearlySingularMatrix');
    
    % Determine indices for k-means cross validation
    KFoldIndices = crossvalind('KFold', Index, numKFolds);
    numPerFold = arrayfun(@(x) nnz(KFoldIndices==x), 1:numKFolds);
    N = max(numPerFold);
    
    % Compute MLR
    currentWeights = zeros(numROIs+1, numStims-1, numKFolds);
    predictions = nan(N,numKFolds,numS);
    for k = 1:numKFolds
        
        % Compute logistic regression
        [currentWeights(:,:,k),~,stats] = mnrfit(Data(:,KFoldIndices~=k)', StimIndex(KFoldIndices~=k,:));
        
        % Validate & generate confusion matrix
        pihat = mnrval(currentWeights(:,:,k), Data(:,KFoldIndices==k)', stats);
        [~,est] = max(pihat,[],2); % determine best guess for each test trial
        predictions(:,k,:) = [est;nan(N-numPerFold(k),1)];
        
        if verbose; parfor_progress; end
    end
    
    pred = zeros(numTrials,numS); % reorder predictions to match StimIndex
    for k = 1:numKFolds
        pred(KFoldIndices==k,:) = predictions(1:numPerFold(k),k,:);
    end
    confusionMatrix(:,:,n) = confusionmat(pred,StimIndex); % compute confusion matrix
    Weights(:,:,n) = mean(currentWeights,3);               % take mean of weights over k-folds
    
end
if verbose; parfor_progress(0); end

% Take mean over repeats
confusionMatrix = mean(confusionMatrix,3);
Weights = mean(Weights,3);


%% Compute stats
stats = struct();
normCM = bsxfun(@rdivide, confusionMatrix, sum(confusionMatrix,1));

% Compute percent correct
stats.PCC = trace(normCM)/numStims;

% Compute selectivity (Elie & Theunissen 2015)
d = diag(normCM);
[Max, maxS] = max(d);
stats.Sel = log2(Max*(numStims-1))/sum(d(setdiff(1:numStims, maxS))); % equal to log2(Sel(rindex)/mean(d(~temp2)));

