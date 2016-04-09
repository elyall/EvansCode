function [dCoM, dSel] = permutationTest(ROIs, varargin)

numPerms = 1000;
distBetween = 1; % CoM

ROIindex = [1,inf];
TrialIndex = [1,inf];

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
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case 'perc'
                perc = varargin{index+1};
                index = index + 2;
            case 'DistBtwn'
                distBetween = varargin{index+1};
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

if ~exist('ROIs', 'var')
    [ROIs, p] = uigetfile({'*.rois'}, 'Select 2 corresponding ROI files', directory,'MultiSelect','on');
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p, ROIs);
end
if numel(ROIs) ~= 2
    error('Requires 2 ROI datasets: Full Pad data and its corresponding Single Whisker data');
end


%% Load data
if iscellstr(ROIs)
    ROIFiles = ROIs;
    for findex = 1:2
        load(ROIFiles, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata; clear ROIdata;
    end
end

% Compute mean stim
for findex = 1:2
    if ~isfield(ROIs{findex}.rois, 'stimMean')
        ROIs{findex} = computeTrialMean(ROIs{findex});
    end
end


%% Determine data to analyze

% Determine ROIs
if isrow(ROIindex)
    ROIindex = ROIindex';
end
if size(ROIindex,2) == 1
    ROIindex = repmat(ROIindex,1,2);
end
if ROIindex(end,1) == inf
    ROIindex = cat(1, ROIindex(1:end-1,:), repmat((ROIindex(1:end-1,1)+1:numel(ROIs{1}.rois))',1,2));
end
numROIs = size(ROIindex,1);

% Determine trials
if ~iscell(TrialIndex)
    TrialIndex = {TrialIndex};
end
if numel(TrialIndex) == 1
    TrialIndex = repmat(TrialIndex,2,1);
end
for findex = 1:2
    if TrialIndex{findex}(end) == inf
        TrialIndex{findex} = cat(2, TrialIndex{findex}(1:end-1), TrialIndex{findex}(1:end-1)+1:numel(ROIs{findex}.DataInfo.StimID));
    end
end
numTrialsPerFile = cellfun(@numel,TrialIndex);


%% Determine permutations

% Determine trials for each stimulus
StimIDs = unique([ROIs{1}.DataInfo.StimID(ismember(ROIs{1}.DataInfo.TrialIndex,TrialIndex{1}));ROIs{2}.DataInfo.StimID(ismember(ROIs{2}.DataInfo.TrialIndex,TrialIndex{2}))]);
numStim = numel(StimIDs);
Trials = cell(numStim,2);
for findex = 1:2
    for sindex = 1:numStim
        Trials{sindex,findex} = find(ROIs{findex}.DataInfo.StimID==StimIDs(sindex));
        Trials{sindex,findex}(~ismember(ROIs{findex}.DataInfo.TrialIndex(Trials{sindex}),TrialIndex{findex})) = [];
    end
end
numTrialsPerStim = cellfun(@numel,Trials);

% Offset and combine Trial indices
Indices = cell(numStim,1);
for sindex = 1:numStim
    Indices{sindex} = [Trials{sindex,1};Trials{sindex,2}+numTrialsPerFile(1)];
end
numTrialsPerStim = cumsum(numTrialsPerStim,2);


%% Compute new data

% Create data matrix of trial means
Data = [ROIs{1}.rois(ROIindex(:,1)).stimMean];
Data = cat(1,Data,[ROIs{2}.rois(ROIindex(:,2)).stimMean]);

% Initialize output
dCoM = nan(numROIs,numPerms);
dSel = nan(numROIs,numPerms);

% Cycle through permutations
parfor_progress(numPerms);
for pindex = 1:numPerms
    
    % Determine permutation
    
    
    % Generate tuning curves (ignore control position)
    Curves = nan(numROIs,numStim-1,2);
    for sindex = 2:numStim
        TI = Indices{sindex}(randperm(numTrialsPerStim(sindex,2)));                     % permute all trials for current stimulus
        Curves(:,sindex-1,1) = mean(Data(TI(1:numTrialsPerStim(sindex,1)),:),1);        % average selected trial means for each ROI
        Curves(:,sindex-1,2) = mean(Data(TI(numTrialsPerStim(sindex,1)+1:end),:),1);    % average remaining trial means for second set of curves
    end
    
    % Compute metrics
    dCoM(:,pindex) = diff([computeCenterOfMass(Curves(:,:,1),1,distBetween),computeCenterOfMass(Curves(:,:,2),1,distBetween)],[],2);
    dSel(:,pindex) = diff([computeVectorSelectivity(Curves,[],1),computeVectorSelectivity(Curves,[],1)],[],2);
    
    parfor_progress; % update status
end
parfor_progress(0);


