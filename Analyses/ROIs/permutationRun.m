function [CoM,p] = permutationRun(ROIs, varargin)

numPerms = 10000;
StimIDs = [];       % vector of indices of stimuli to analyze
ROIindex = [];      % matrix of dimensions numROIs by numFiles
TrialIndex = {};    % cell of length numFiles, where each element is indices of trials
outliers = {};      % cell of length numFiles, where each element is totalTrials x totalROIs

% CoM calculation
positions = [];
distBetween = 1;

saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numPerms'
                numPerms = varargin{index+1};
                index = index + 2;
            case 'StimIDs'
                StimIDs = varargin{index+1};
                index = index + 2;
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'TrialIndex', 'Trials', 'trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'outliers'
                outliers = varargin{index+1};
                index = index + 2;
            case 'positions'
                positions = varargin{index+1};
                index = index + 2;
            case {'DistBtwn','distBetween'}
                distBetween = varargin{index+1};
                index = index + 2;
            case {'save','Save'}
                saveOut = true;
                index = index + 1;
            case 'saveFile'
                saveFile = varargin{index+1};
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
if isstruct(ROIs)
    ROIs = {ROIs};
end


%% Load data
numFiles = numel(ROIs);
if iscellstr(ROIs)
    ROIFiles = ROIs;
    for findex = 1:numFiles
        load(ROIFiles{findex}, 'ROIdata', 'outliers', '-mat');
        ROIs{findex} = ROIdata; clear ROIdata;
    end
    if saveOut && isempty(saveFile)
        saveFile = strcat(strtok(ROIFiles{1},'.'),'.perm');
    end
end

% Compute mean stim
for findex = 1:numFiles
    if ~isfield(ROIs{findex}.rois, 'stimMean')
        ROIs{findex} = computeTrialMean(ROIs{findex});
    end
end


%% Determine ROIs
totalROIs = cellfun(@(x) numel(x.rois), ROIs);
if isempty(ROIindex)
    ROIindex = cell(1,numFiles);
    for findex = 1:numFiles
        ROIindex{findex} = (1:totalROIs(findex))';
    end
    try
        ROIindex = cat(2,ROIindex{:});
    catch
        error('Files need to have the same number of ROIs, or ROIindex needs to be specified');
    end
else % ROIindex specified
    if isrow(ROIindex)
        ROIindex = ROIindex';
    end
    if size(ROIindex,2) ~= numFiles
        ROIindex = repmat(ROIindex,1,2); % assumes index is same for both files
    end
end
numROIs = size(ROIindex,1);


%% Determine trials
totalTrials = cellfun(@(x) numel(x.DataInfo.TrialIndex), ROIs);
if isempty(TrialIndex) % use all trials
    TrialIndex = cell(1,numFiles);
    for findex = 1:numFiles
        TrialIndex{findex} = true(totalTrials(findex),1);
    end
else
    if ~iscell(TrialIndex)
        TrialIndex = {TrialIndex};
    end
    if numel(TrialIndex) ~= numFiles
        TrialIndex = repmat(TrialIndex,1,2); % assumes index is same for both files
    end
    for findex = 1:numFiles
        TrialIndex{findex} = ismember(ROIs{findex}.DataInfo.TrialIndex',TrialIndex{findex}); % convert to logical vector
    end
end
numTrials = cellfun(@nnz, TrialIndex);


%% Determine outliers
if isempty(outliers)
    outliers = cell(1,numFiles);
    for findex = 1:numFiles
        outliers{findex} = false(totalTrials(findex),totalROIs(findex));
    end
elseif isequal(outliers, true) % compute outliers
    outliers = cell(1,numFiles);
    for findex = 1:numFiles
        outliers{findex} = false(totalTrials(findex),totalROIs(findex));
        for rindex = ROIindex(:,findex)';
            outliers{findex}(TrialIndex{findex},rindex) = determineOutliers(ROIs{findex}.rois(rindex).stimMean(TrialIndex{findex}),'GroupID',ROIs{findex}.DataInfo.StimID(TrialIndex{findex}),'type','medianRule');
        end
    end
else % outliers input
    if ~iscell(outliers)
        outliers = {outliers};
    end
    if numel(outliers) ~= numFiles
        outliers = repmat(outliers,1,2);
    end
end


%% Determine stimuli
if isempty(StimIDs)
    StimIDs = cellfun(@(x) unique(x.DataInfo.StimID),ROIs,'UniformOutput',false);
    StimIDs = unique(cat(1,StimIDs{:}));
    StimIDs(1) = []; % remove control stimulus
end
numStim = numel(StimIDs);

% Determine trials for each stimulus
StimIndex = cell(1,numFiles);
for findex = 1:numFiles
    StimIndex{findex} = false(totalTrials(findex),numStim); % empty dimension to make bsxfun step easier
    for sindex = 1:numStim
        StimIndex{findex}(:,sindex) = ROIs{findex}.DataInfo.StimID==StimIDs(sindex);
    end
end


%% Gather data to analyze

% Determine trials to use for each ROI
Trials = cell(1,numFiles);
for findex = 1:numFiles
    Trials{findex} = bsxfun(@and, TrialIndex{findex}, ~outliers{findex}(:,ROIindex(:,findex)));
    Trials{findex} = bsxfun(@and, Trials{findex}, permute(StimIndex{findex},[1,3,2]));
end
if numFiles == 2
    numTrialsPerStim = cellfun(@sum,Trials,'UniformOutput',false);
    numTrialsPerStim = squeeze(cat(4, numTrialsPerStim{:}));
    numTrialsPerStim = cumsum(numTrialsPerStim,3);
elseif numFiles == 1
    midwaypoint = floor(totalTrials/2);
    numTrialsPerStim = squeeze(sum(Trials{1}(1:midwaypoint,:,:)));
    numTrialsPerStim = cat(3, numTrialsPerStim, squeeze(sum(Trials{1})));
end
Trials = cat(1,Trials{:});

% Gather Data: matrix of trial means
StimMeans = [];
for findex = 1:numFiles
    StimMeans = cat(1,StimMeans,[ROIs{findex}.rois(ROIindex(:,findex)).stimMean]);
end
num = max(max(numTrialsPerStim(:,:,2)));
Data = nan(numStim,numROIs,num);
for rindex = 1:numROIs
    for sindex = 1:numStim
        Data(sindex,rindex,:) = [StimMeans(Trials(:,rindex,sindex),rindex);nan(num-numTrialsPerStim(rindex,sindex,2),1)];
    end
end


%% Compute metrics off of permutations

% Initialize output
dCoM = nan(numROIs,numPerms+1);
dSel = nan(numROIs,numPerms+1);

% Cycle through permutations
fprintf('Computing %d permutations for %d ROIs...\n',numPerms,numROIs);
pfH = parfor_progress(numPerms+1);
parfor pindex = 1:numPerms+1
    
    % Generate tuning curves (ignore control position)
    Curves = nan(numROIs,numStim,2);
    for sindex = 1:numStim
        for rindex = 1:numROIs
            if pindex ~= 1
                TI = randperm(numTrialsPerStim(rindex,sindex,2));
            else
                TI = 1:numTrialsPerStim(rindex,sindex,2);
            end
            Curves(rindex,sindex,1) = mean(Data(sindex,rindex,TI(1:numTrialsPerStim(rindex,sindex,1))));
            Curves(rindex,sindex,2) = mean(Data(sindex,rindex,TI(numTrialsPerStim(rindex,sindex,1)+1:numTrialsPerStim(rindex,sindex,2))));
        end
    end
    
    % Compute center of mass change
    dCoM(:,pindex) = diff([computeCenterOfMass(Curves(:,:,1),positions,distBetween),computeCenterOfMass(Curves(:,:,2),positions,distBetween)],[],2);
    
    % Compute selectivity change
    Min = min(min(Curves,[],3),[],2);
    dSel(:,pindex) = diff([computeVectorSelectivity(Curves(:,:,1),Min,1),computeVectorSelectivity(Curves(:,:,2),Min,1)],[],2);
    
    parfor_progress(pfH); % update status
end
parfor_progress(pfH,0);

Actual = [dCoM(:,1),dSel(:,1)];
dCoM(:,1) = [];
dSel(:,1) = [];


%% Generate p-values
p = nan(numROIs,2);
p(:,1) = sum(bsxfun(@gt,abs(dCoM),abs(Actual(:,1))),2)/numPerms;
p(:,2) = sum(bsxfun(@gt,abs(dSel),abs(Actual(:,2))),2)/numPerms;


%% Save outputs
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile,'file')
        save(saveFile,'dCoM','dSel','p','Actual','ROIindex','-v7.3');
    else
        save(saveFile,'dCoM','dSel','p','Actual','ROIindex','-append');
    end
    fprintf('Saved permutation results to: %s\n',saveFile);
end


