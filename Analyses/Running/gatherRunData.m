function [RunData, TrialIndex, StimID] = gatherRunData(AnalysisInfo,frames,varargin)

TrialIndex = [1 inf];
type = 'mat'; % 'mat' or 'cell'

directory = cd;

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
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if (~exist('AnalysisInfo','var') || isempty(AnalysisInfo)) && (~exist('frames', 'var') || isempty(frames))
    [AnalysisInfo, p] = uigetfile({'*.exp'},'Choose Experiment file',directory);
    if isnumeric(AnalysisInfo)
        return
    end
    AnalysisInfo = fullfile(p,AnalysisInfo);
end


%% Load data
if ischar(AnalysisInfo)
    load(AnalysisInfo,'AnalysisInfo','frames','-mat');
end


%% Determine data to pull out
if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1),TrialIndex(1:end-1)+1:size(AnalysisInfo,1)];
end
numTrials = numel(TrialIndex);


%% Pull out running data
RunData = cell(numTrials,1);
for tindex = 1:numTrials
    if exist('frames','var') && ~isempty(frames)
        RunData{tindex} = frames.RunningSpeed(AnalysisInfo.ExpStimFrames(TrialIndex(tindex),1):AnalysisInfo.ExpStimFrames(TrialIndex(tindex),2));
    else
        RunData{tindex} = AnalysisInfo.meanRunningSpeed{TrialIndex(tindex)}(AnalysisInfo.TrialStimFrames(TrialIndex(tindex),1):AnalysisInfo.TrialStimFrames(TrialIndex(tindex),2));
    end
end
StimID = AnalysisInfo.StimID(TrialIndex);


%% Format output
if strcmp(type,'mat') % Convert Data to a matrix
    D = cellfun(@numel,RunData);
    temp = nan(numTrials,max(D));
    for tindex = 1:numTrials
        temp(tindex,1:D(tindex)) = RunData{tindex};
    end
    RunData = temp;
end