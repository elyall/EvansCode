function [Probabilities, numTrials, StimIDs] = NeuroMetricFunction(ROIdata, ROIid, FrameIndex, varargin)


type = 'dFoF'; % 'dFoF' or 'spikes'
thresh = .2;


%% Parse input arguments
if ~exist('ROIdata','var') || isempty(ROIdata)
    directory = CanalSettings('DataDirectory');
    [ROIdata, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if ~exist('ROIid', 'var') || isempty(ROIid)
    ROIid = [1 inf];
elseif ischar(ROIid) && strcmp(ROIid, 'all')
    ROIid = [1 inf];
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case {'Threshold', 'threshold', 'thresh'}
                thresh = varargin{index+1};
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


%% Load in data
if ischar(ROIdata)
    load(ROIdata, 'ROIdata');
end
StimIDs = unique(ROIdata.DataInfo.StimID);
numStims = numel(StimIDs);


%% Determine frames to analyze
if ~exist('FrameIndex', 'var') || isempty(FrameIndex)
    FrameIndex = [ROIdata.DataInfo.numFramesBefore+1, ROIdata.DataInfo.numFramesBefore + mode(ROIdata.DataInfo.numStimFrames)];
end


%% Determine ROIs to filter
if ROIid(end) == inf
    ROIid = [ROIid(1:end-1), ROIid(1:end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIid);


%% Calculate probability of going above threshold for each stimulus

Probabilities = nan(numROIs, numStims);
numTrials = zeros(1,numStims);

for sindex = 1:numStims
    TrialIndex = ROIdata.DataInfo.StimID == StimIDs(sindex);
    numTrials(sindex) = sum(TrialIndex);
    for rindex = 1:numROIs
        switch type
            case 'dFoF'
                data = ROIdata.rois(rindex).dFoF(TrialIndex, FrameIndex(1):FrameIndex(2));
            case 'spikes'
                data = ROIdata.rois(rindex).spike(TrialIndex, FrameIndex(1):FrameIndex(2));
        end
        index = any(data>=thresh, 2);
        Probabilities(rindex, sindex) = sum(index)/numTrials(sindex);
    end
end

