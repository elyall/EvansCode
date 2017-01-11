function CoM = computeTCOverTime(ROIdata, duration, startOfWindow, varargin)

saveOut = false;
saveFile = '';

ROIindex = [1 inf];
TrialIndex = [1 inf];
frameRate = 15.46;
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
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
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

if ~exist('duration','var') || isempty(duration)
    duration = .5; % seconds
end

if ~exist('startOfWindow','var') || isempty(startOfWindow)
    startOfWindow = 0:.1:1; % seconds
end
numWindows = numel(startOfWindow);

fprintf('Calculating CoM with sliding window over stimulus...\n');


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
numTrials = size(ROIdata.rois(1).dFoF,1);


%% Determine data to analyze
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIindex);

if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:numTrials];
end

% Determine column indices of windows
numFrames = round(frameRate*duration);
FrameIndex = ROIdata.DataInfo.numFramesBefore+1:ROIdata.DataInfo.numFramesBefore+numFrames; % column indices for first window within stimulus
maxStimFrames = max(ROIdata.DataInfo.numStimFrames);
stimFrameTimeStamps = 0:1/frameRate:(maxStimFrames-1)/frameRate;            % time of start of each frame 
[~,Offset] = min(abs(bsxfun(@minus, startOfWindow,stimFrameTimeStamps')));  % determine first frame that overlaps with start of window by >=50% of that frame
FrameIndex = bsxfun(@plus, FrameIndex,(Offset-1)');                         % column indices for each window


%% Calculate center of mass via sliding window
CoM = nan(numROIs, numWindows);
parfor windex = 1:numWindows
    temp = computeTrialMean(ROIdata, 'ROIindex', ROIindex, 'FrameIndex', FrameIndex(windex,:));
    [~, Curves, ~] = computeTuningCurve(temp, ROIindex, TrialIndex, 'ControlID', ControlID, 'StimIDs', StimIDs);
    CoM(:,windex) = computeCenterOfMass(Curves, positions, distBetween);
end


%% Save output
if saveOut
    if exist(saveFile,'file')
        save(saveFile,'CoM','-append');
    else
        save(saveFile,'CoM','-v7.3');
    end
    fprintf('Saved sliding window analysis to: %s\n',saveFile);
end


%% Display results
if verbose
    figure;
    errorbar(startOfWindow,mean(CoM),std(CoM)/sqrt(numROIs));
    xlabel(sprintf('Start of %dms window',duration*1000));
    ylabel('Center of Mass');
end

