function CoM = computeTCOverTime(ROIdata, duration, shift, varargin)

fn = '~/Documents/AdesnikLab/Data/2481_179_000_final.rois';

saveOut = false;
saveFile = '';

ROIindex = [1 inf];
TrialIndex = [1 inf];
FrameIndex = [];

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
            case 'FrameIndex'
                FrameIndex = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
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

if ~exist('shift','var') || isempty(shift)
    shift = .1; % seconds
end

fprintf('Calculating CoM with sliding window over stimulus...');


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

if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:numTrials];
end

if isempty(FrameIndex)
    FrameIndex = cell(numTrials,1);
    for tindex = 1:numTrials
%         FrameIndex{tindex} = ROIdata.DataInfo.numFramesBefore+1:ROIdata.DataInfo.numFramesBefore+ROIdata.DataInfo.numStimFrames(tindex);
        F = ROIdata.DataInfo.numFramesBefore+ROIdata.DataInfo.numStimFrames(tindex);
        FrameIndex{tindex} = F-7:F;
    end
elseif ~iscell(FrameIndex)
    FrameIndex = {FrameIndex};
end
if numel(FrameIndex) == 1
    FrameIndex = repmat(FrameIndex,numTrials,1);
end


%% Calculate via sliding window
ControlID = [];
StimIDs = [];
positions = [];
distBetween = [];

FrameIndex = ROIdata.DataInfo.numFramesBefore+1:ROIdata.DataInfo.numFramesBefore+round(frameRate*duration); % initial window of frames
numFrames = mode(ROIdata.DataInfo.numStimFrames(tindex));
Offset = round(0:1/frameRate:numFrames/frameRate);
[~,Offset] = min(abs(bsxfun(@minus, Offset, 0:.1:numFrames)));

CoM = nan(numROIs, N);
for tindex = 1:N
    currentFrames = FrameIndex+Offset(tindex);
    ROIdata = computeTrialMean(ROIdata, 'ROIindex', ROIindex, 'FrameIndex', currentFrames);
    [~, Curves, ~] = computeTuningCurve(ROIdata, ROIindex, TrialIndex, 'ControlID', ControlID, 'StimIDs', StimIDs);
    CoM(:,N) = computeCenterOfMass(Curves, positions, distBetween);
end

