function [TrialData,numFramesBefore,numStimFrames] = trialOrganize(Data,AnalysisInfo,depthID,varargin)
% Data is totalFrames x numROIs
% convert vector of all data to numTrials x numFrames x numROIs

numFramesBefore = 16;
numFramesAfter = 31;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numFramesBefore'
                numFramesBefore = varargin{index+1};
                index = index + 2;
            case 'numFramesAfter'
                numFramesAfter = varargin{index+1};
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

%%
[totalFrames,numROIs] = size(Data);
numTrials = size(AnalysisInfo,1);
numFrames = numFramesBefore+numFramesAfter+1;

StimFrames = AnalysisInfo.ExpStimFrames;
FirstFrame = arrayfun(@(x) find(depthID>=StimFrames(x,1),1,'first'), 1:numTrials);
LastFrame = arrayfun(@(x) find(depthID<StimFrames(x,2),1,'last'), 1:numTrials);
numStimFrames = LastFrame - FirstFrame + 1;

FrameIndex = [max(FirstFrame-numFramesBefore,1)',min(FirstFrame+numFramesAfter,totalFrames)'];
TrialData = nan(numFrames,numTrials,numROIs);
for r = 1:numROIs
    for t = 1:numTrials
        current = Data(FrameIndex(t,1):FrameIndex(t,2),r);
        if numel(current) == numFrames
            TrialData(:,t,r) = current;
        else
            temp = nan(numFrames-numel(current),1);
            if t==1
                TrialData(:,t,r) = [temp;current];
            else
                TrialData(:,t,r) = [current;temp];
            end
        end
    end
end

fprintf('\tComplete\n');