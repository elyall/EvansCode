function [CoM,Max,p_tuned,numTrialsPerStim,StimIDs] = computeTCOverMeanSpeed(ROIdata, TrialRunSpeed, width, startOfWindow, varargin)

saveOut = false;
saveFile = '';

ROIindex = [1 inf];
TrialIndex = [1 inf];
verbose = false;

% Tuning curve parameters
ControlID = 0;
StimIDs = [];

% Center of mass parameters
positions = 2;
distBetween = 1;

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
            case 'ControlID'
                ControlID = varargin{index+1};
                index = index + 2;
            case 'StimIDs'
                StimIDs = varargin{index+1};
                index = index + 2;
            case 'positions'
                positions = varargin{index+1};
                index = index + 2;
            case 'distBetween'
                distBetween = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
                index = index + 1;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('ROIdata','var') || isempty(ROIdata)
    [ROIdata, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if ~exist('TrialRunSpeed','var') || isempty(TrialRunSpeed)
    TrialRunSpeed = closestFile(ROIdata.filename,'.exp');
    TrialRunSpeed = TrialRunSpeed{1};
end
if isempty(TrialRunSpeed)
    error('need run speed data');
end

if ~exist('width','var') || isempty(width)
    width = 100; % deg/s or cm/s (depends on TrialRunSpeed input)
end

if ~exist('startOfWindow','var') || isempty(startOfWindow)
    startOfWindow = 0:25:300; % deg/s or cm/s (depends on TrialRunSpeed input)
end
numWindows = numel(startOfWindow);

fprintf('Calculating CoM with sliding window over mean speed...\n');


%% Load ROI data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end

% Compute dF/F
if ~isfield(ROIdata.rois, 'dFoF')
    ROIdata = computeDFoF(ROIdata);
end

% Compute trial mean
if ~isfield(ROIdata.rois, 'stimMean')
    ROIdata = computeTrialMean(ROIdata);
end


%% Determine data to analyze
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIindex);

if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:ROIdata.DataInfo.TrialIndex(end)];
end
TrialIndex = ismember(ROIdata.DataInfo.TrialIndex, TrialIndex);
numTrials = numel(TrialIndex);


%% Gather run speed data
if ischar(TrialRunSpeed)
    TrialRunSpeed = gatherRunData(TrialRunSpeed);
end
SpeedMean = nanmean(TrialRunSpeed(TrialIndex,:),2);
% SpeedStd = nanstd(TrialRunSpeed(TrialIndex),2);

% Detrmine trials in each window
Index = false(numTrials,numWindows);
for windex = 1:numWindows
    Index(:,windex) = SpeedMean>=startOfWindow(windex) & SpeedMean<=startOfWindow(windex)+width;
end

% Determine number of trials per stimulus
if isempty(StimIDs)
    StimIDs = unique(ROIdata.DataInfo.StimID(TrialIndex))';
end
numTrialsPerStim = nan(numel(StimIDs),numWindows);
for windex = 1:numWindows
    for sindex = 1:numel(StimIDs)
        numTrialsPerStim(sindex,windex) = nnz(ROIdata.DataInfo.StimID(TrialIndex' & Index(:,windex))==StimIDs(sindex));
    end
end

% Ignore windows without a single trial for any given stimulus
bad = any(numTrialsPerStim==0,1);


%% Calculate center of mass via sliding window
CoM = nan(numROIs, numWindows);
p_tuned = nan(numROIs, numWindows);
Max = nan(numROIs, numWindows);
TrialIndex = ROIdata.DataInfo.TrialIndex(TrialIndex);
parfor windex = find(~bad)
    [~, Curves, ~, p_tuned(:,windex)] = computeTuningCurve(ROIdata, ROIindex, TrialIndex(Index(:,windex)), 'ControlID', ControlID, 'StimIDs', StimIDs);
    Max(:,windex) = max(Curves(:,2:end),[],2);
    CoM(:,windex) = computeCenterOfMass(Curves, positions, distBetween);
end


%% Save output
if saveOut
    if exist(saveFile,'file')
        save(saveFile,'CoM','Max','p_tuned','numTrialsPerStim','StimIDs','-append');
    else
        save(saveFile,'CoM','Max','p_tuned','numTrialsPerStim','StimIDs','-v7.3');
    end
    fprintf('Saved sliding window analysis to: %s\n',saveFile);
end


%% Display results
if verbose
    
    % Compute mean and standard error for significantly tuned and driven neurons
    Mean = nan(1,numWindows);
    StdErrMean = nan(1,numWindows);
    for windex = 1:numWindows
        index = p_tuned(:,windex)<0.05 & Max(:,windex)>.2;
        Mean(windex) = mean(CoM(index,windex));
        StdErrMean(windex) = std(CoM(index,windex))/sqrt(nnz(index));
    end
    
    % Display curve
    figure;
    errorbar(startOfWindow,mean(CoM),std(CoM)/sqrt(numROIs));
    xlabel(sprintf('Start of %dms window',width*1000));
    ylabel('Center of Mass');
end

