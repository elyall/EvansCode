function [p,Corr,CoM,TrialIndex,outliers] = permutationRun(ROIs, targetTrials, varargin)

numPerms = 10000;
StimIDs = [];       % vector of indices of stimuli to analyze
ROIindex = [];     
TrialIndex = [];    % indices of trials to analyze
RunSpeed = [];      % mean run speed for each trial
outliers = [];      % logical matrix specifying which trials are outliers for a given ROI (totalTrials x totalROIs)

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
    [ROIs, p] = uigetfile({'*.rois'}, 'Select 2 corresponding ROI files', directory);
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p, ROIs);
end


%% Load data
if ischar(ROIs)
    ROIFile = ROIs;
    for findex = 1:numFiles
        load(ROIFile{findex}, 'ROIdata', 'outliers', '-mat');
        ROIs = ROIdata; clear ROIdata;
    end
    if saveOut && isempty(saveFile)
        saveFile = strcat(strtok(ROIFile{1},'.'),'.permrun');
    end
end

% Compute mean stim
if ~isfield(ROIs.rois, 'stimMean')
    ROIs = computeTrialMean(ROIs{findex});
end


%% Determine ROIs
totalROIs = numel(ROIs.rois);
if isempty(ROIindex)
    ROIindex = 1:totalROIs;
end
numROIs = numel(ROIindex);


%% Determine trials
totalTrials = numel(ROIs.DataInfo.TrialIndex);
if isempty(TrialIndex) % use all trials
    TrialIndex = true(totalTrials,1);
else
    TrialIndex = ismember(ROIs.DataInfo.TrialIndex',TrialIndex); % convert to logical vector
end
numTrials = nnz(TrialIndex);
targetTrials = ismember(ROIs.DataInfo.TrialIndex',targetTrials);


%% Determine outliers
if isempty(outliers)
    outliers = false(totalTrials,numROIs);
elseif isequal(outliers, true) % compute outliers
    outliers = false(totalTrials,totalROIs);
    for rindex = ROIindex'
        outliers(TrialIndex,rindex) = determineOutliers(ROIs.rois(rindex).stimMean(TrialIndex),'GroupID',ROIs.DataInfo.StimID(TrialIndex),'type','medianRule');
    end
end


%% Determine stimuli
if isempty(StimIDs)
    StimIDs = unique(ROIs.DataInfo.StimID);
    StimIDs(1) = []; % remove control stimulus
end
numStim = numel(StimIDs);

% Determine trials for each stimulus
StimIndex = bsxfun(@eq, ROIs.DataInfo.StimID, StimIDs');

% Determine trials available for each stimulus
Trials = bsxfun(@and, TrialIndex, StimIndex);
totalTrialsPerStim = sum(Trials);


%% Gather data to analyze

% Determine number of trials to pull for each ROI for each stimulus
targetTrials = bsxfun(@and, TargetTrials, StimIndex);
numTrialsPerStim = sum(targetTrials);

% Gather Data: matrix of trial means
StimMeans = repmat(permute([ROIs.rois(ROIindex).stimMean],[1,3,2]),1,numStim,1);
StimMeans(~repmat(StimIndex,1,numROIs,1)) = NaN;   % remove data for other stimuli
StimMeans(repmat(permute(outliers(:,ROIindex),[1,3,2]),1,numStim,1)) = NaN;          % remove outliers


%% Compute metrics off of permutations

% Initialize output
CoM = nan(numROIs,numPerms+1);
TrialIndex = false(totalTrials,numPerms+1);

% Cycle through permutations computing map
fprintf('Computing CoM for %d ROIs (%d permutations)...\n',numROIs,numPerms);
pfH = parfor_progress(numPerms+1);
parfor pindex = 1:numPerms+1
    
    % Determine trials to sample
    if pindex == 1
        currentTrials = repmat(targetTrials,1,1,numROIs);
    else
        currentTrials = Trials; % FIX
    end
    TrialIndex(:,pindex) = any(currentTrials,2);
    
    % Generate tuning curves (ignore control position)
    Curves = nanmean(StimMeans(repmat(currentTrials,1,1,numROIs)));
    
    % Compute center of mass 
    CoM(:,pindex) = computeCenterOfMass(Curves,positions,distBetween);
    
    
    parfor_progress(pfH); % update status
end
parfor_progress(pfH,0);


%% Determine runspeed mean and variance for each permutation
if ~isempty(RunSpeed)
    RunSpeed = repmat(RunSpeed,1,numPerms+1);
    Speed.mean = mean(RunSpeed(TrialIndex));
    Speed.std = std(RunSpeed(TrialIndex));
else
    Speed = [];
end


%% Compute correlations of map
Corr = nan(1,numPerms+1);
Loc = [ROIs.rois(ROIindex).centroid];

fprintf('Computing map correlation for %d ROIs (%d permutations)...\n',numROIs,numPerms);
pfH = parfor_progress(numPerms+1);
parfor pindex = 1:numPerms+1
    
    % Compute map correlation
    Corr(pindex) = computeCorr(CoM(:,pindex),Loc); %FIX
    
    parfor_progress(pfH); % update status
end
parfor_progress(pfH,0);

% Generate p-value
temp = mean(Corr);
Corr = Corr - temp;
p = sum(abs(Corr(2:end))>abs(Corr(1)))/numPerms;
Corr = Corr + temp;


%% Save outputs
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile,'file')
        save(saveFile,'dCoM','dSel','p','Actual','ROIindex','-v7.3');
    else
        save(saveFile,'dCoM','dSel','p','Actual','ROIindex','-append');
    end
    fprintf('Saved permutation results to: %s\n',saveFile);
end


