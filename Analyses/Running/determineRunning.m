function [TrialIndex, RunIndex] = determineRunning(AnalysisInfo, frames, thresh, varargin)

directory = cd;
TrialIndex = [1 inf];


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'TrialIndex'
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case 'comparison'
                comparison = varargin{index+1};
                index = index + 2;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if (~exist('AnalysisInfo', 'var') && ~exist('frames', 'var')) || (isempty(AnalysisInfo) && isempty(frames))
    [ExperimentFile, p] = uigetfile('Select Experiment file:', directory);
    if ~ExperimentFile
        return
    end
    ExperimentFile = fullfile(p, ExperimentFile);
elseif ischar(AnalysisInfo)
    ExperimentFile = AnalysisInfo;
    clear AnalysisInfo
end

if ~exist('thresh', 'var') || isempty(thresh)
    thresh = 100; % degrees per second
end

%% Load in variables
if ~exist('AnalysisInfo', 'var') || isempty(AnalysisInfo)
    load(ExperimentFile, 'AnalysisInfo', '-mat');
end
if ~exist('frames', 'var') || isempty(frames)
    load(ExperimentFile, 'frames', '-mat');
end


%% Determine trials to analyze
if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:size(AnalysisInfo, 1)];
end
numTrials = numel(TrialIndex);


%% Determine in which trials the mouse is running

% initialize output
RunIndex = false(numTrials, 1);
    
for tindex = TrialIndex
    if mean(frames.RunningSpeed(AnalysisInfo.ExpStimFrames(tindex,1):AnalysisInfo.ExpStimFrames(tindex,2))) > thresh
        RunIndex(tindex) = true;
    end
end    

% Remove non-running trials from index
TrialIndex(~RunIndex) = [];

