function [Weights, confusionMatrix, stats] = computeMLR(Data, StimIndex, varargin)
% Data is numROIs by numTrials

numKFolds = 5;
% ControlID = [];

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
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Data', 'var')
    [Data, p] = uigetfile({'*.rois;*.mat'},'Select ROI file',directory,'MultiSelect','on');
    if ~Data
        return
    else
        Data = fullfile(p, Data);
    end
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
[numROIs,numTrials] = size(Data);

% Determine stimuli indices
[StimIDs,~,StimIndex] = unique(StimIndex);
numStims = numel(StimIDs);


%% Perform analysis

% Determine indices for k-means cross validation
KFoldIndices = crossvalind('KFold', StimIndex, numKFolds);

% Initialize output
Weights = zeros(2, numStims-1, numROIs);
confusionMatrix = zeros(numStims, numStims, numROIs);

% Compute MLR
parfor_progress(numROIs*numKFolds);
% warning('off', 'stats:mnrfit:IterOrEvalLimit');
% warning('off', 'MATLAB:nearlySingularMatrix');
parfor rindex = 1:numROIs
    predictions = zeros(numTrials,1);
    for k = 1:numKFolds
        
        % Compute logistic regression
        [currentWeights,~,Statistics] = mnrfit(Data(rindex,KFoldIndices~=k), StimIndex(KFoldIndices~=k));
        Weights(:,:,rindex) = Weights(:,:,rindex)*(k-1)/k + currentWeights/k;
        
        % Validate & generate confusion matrix
        pihat = mnrval(currentWeights, Data(rindex,KFoldIndices==k)', Statistics);
        [~,est] = max(pihat,[],2); % determine best guess for each test trial
        predictions(KFoldIndices==k) = est;

        parfor_progress;
    end
    
    confusionMatrix(:,:,rindex) = confusionmat(predictions,StimIndex);
    
end
parfor_progress(0);

Weights = squeeze(Weights(2,:,:));

%% Compute stats

% Compute percent correct & selectivity
PCC = nan(numROIs,1);
Sel = nan(numROIs,1);
parfor rindex = 1:numROIs
    
    % Normalize confusion matrix
    normCM = bsxfun(@rdivide, confusionMatrix(:,:,rindex), sum(confusionMatrix(:,:,rindex),1));
    
    % Compute PCC
    PCC(rindex) = trace(normCM)/numStims;
    
    % Compute Selectivity
    d = diag(normCM);
    [Max, maxS] = max(d);
    Sel(rindex) = log2(Max*(numStims-1)/sum(d(setdiff(1:numStims, maxS))));
    
end

% Create output
stats.PCC = PCC;
stats.Sel = Sel;

