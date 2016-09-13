function [numTrials, StimIDs] = determineNumTrialsPerStim(ROIs,TrialIndex)


%% Load in data
if ~iscell(ROIs)
    ROIs = {ROIs};
end
numFiles = numel(ROIs);
if iscellstr(ROIs)
    ROIFiles = ROIs;
    for findex = 1:numFiles
        load(ROIFiles{findex},'ROIdata','-mat');
        ROIs{findex} = ROIdata; clear ROIdata;
    end
end

if ~exist('TrialIndex','var') || isempty(TrialIndex)
    TrialIndex = cellfun(@(x) x.DataInfo.TrialIndex,ROIs,'UniformOutput',false);
elseif ~iscell(TrialIndex)
    TrialIndex = {TrialIndex};
end
for findex = 1:numFiles
    TrialIndex{findex} = ismember(ROIs{findex}.DataInfo.TrialIndex,TrialIndex{findex});
end


%% Determine stimuli
StimIDs = cellfun(@(x) unique(x.DataInfo.StimID),ROIs,'UniformOutput',false);
numStims = cellfun(@numel,StimIDs);
temp = nan(numFiles,max(numStims));
for findex = 1:numFiles
    temp(findex,1:numStims(findex)) = StimIDs{findex};
end
StimIDs = temp;


%% Determine number of trials
numTrials = nan(size(StimIDs));
for findex = 1:numFiles
    for sindex = 1:numStims(findex)
        numTrials(findex,sindex) = nnz(ROIs{findex}.DataInfo.StimID(TrialIndex{findex})==StimIDs(findex,sindex));
    end
end


