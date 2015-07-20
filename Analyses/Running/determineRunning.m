function TrialIndex = determineRunning(AnalysisInfo, frames, thresh)

directory = cd;

%% Parse input arguments
% index = 1;
% while index<=length(varargin)
%     try
%         switch varargin{index}
%             case 'logical'
%                 outType = 'logical';
%                 index = index + 1;
%             case 'index'
%                 outType = 'index';
%                 index = index + 1;
%         end
%     catch
%         warning('Argument %d not recognized',index);
%         index = index + 1;
%     end
% end

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
numTrials = size(AnalysisInfo, 1);
TrialIndex = false(numTrials, 1);
for tindex = 1:numTrials
    if mean(frames.RunningSpeed(AnalysisInfo.ExpStimFrames(tindex,1):AnalysisInfo.ExpStimFrames(tindex,2))) >= thresh
        TrialIndex(tindex) = true;
    end
end

