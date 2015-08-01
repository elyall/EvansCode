function RunIndex = determineRunning(AnalysisInfo, frames, thresh, varargin)

directory = cd;
type = 'TrialStimMean'; % 'TrialStimMean' or 'ExpStimVar'
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


%% Format running speed
numTrials = size(AnalysisInfo, 1);
numFrames = max(AnalysisInfo.nFrames);
RunningSpeed = nan(numTrials, numFrames);
for tindex = 1:numTrials
    RunningSpeed(tindex, :) = frames.RunningSpeed(AnalysisInfo.ExpFrames(tindex, 1):AnalysisInfo.ExpFrames(tindex, 1) + numFrames -1);
end


%% Determine in which trials the mouse is running

% initialize output
RunIndex = false(numTrials, 1);

switch type
    
    case 'TrialStimMean'
        for tindex = 1:numTrials
            if eval(sprintf('mean(RunningSpeed(tindex, AnalysisInfo.TrialStimFrames(tindex,1):AnalysisInfo.TrialStimFrames(tindex,2))) %s thresh', comparison))
                RunIndex(tindex) = true;
            end
        end
        
    case 'ExpStimThresh'
        StimIDs = unique(AnalysisInfo.StimID);
        figure;
        for sindex = 1:numel(StimIDs)
            currentTrials = find(AnalysisInfo.StimID==StimIDs(sindex));
            stimRunning = RunningSpeed(currentTrials, AnalysisInfo.TrialStimFrames(currentTrials(1), 1):AnalysisInfo.TrialStimFrames(currentTrials(1), 2));
            RunIndex(currentTrials) = mean(stimRunning, 2) >= thresh;
            
            goodTrials = all(stimRunning >= thresh, 2); % all trials for current stimulus in which mouse was always running faster than thresh
            
            subplot(1, 9, sindex); plot(stimRunning(~goodTrials, :)', 'k'); hold on; plot(stimRunning(goodTrials, :)', 'r'); ylim([-200, 700]);
        end
        disp('');
        
    case 'ExpVar'
        
        
end
