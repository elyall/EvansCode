function RunIndex = determineRunning(AnalysisInfo, frames, thresh, varargin)

directory = cd;
type = 'StimPeriods'; % 'StimPeriods' or 'frames'
comparison = '>=';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
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

%% Load in variables
if ~exist('AnalysisInfo', 'var') || isempty(AnalysisInfo)
    load(ExperimentFile, 'AnalysisInfo', '-mat');
end
if ~exist('frames', 'var') || isempty(frames)
    load(ExperimentFile, 'frames', '-mat');
end

%% Determine in which trials the mouse is running
switch type
    case 'StimPeriods'
        numTrials = size(AnalysisInfo, 1);
        RunIndex = false(numTrials, 1);
        for tindex = 1:numTrials
            if eval(sprintf('frames.RunningSpeed(AnalysisInfo.ExpStimFrames(tindex,1):AnalysisInfo.ExpStimFrames(tindex,2)) %s thresh', comparison))
                RunIndex(tindex) = true;
            end
        end
    case 'frames'
        eval(sprintf('RunIndex = frames.RunningSpeed %s thresh;', comparison));
end
