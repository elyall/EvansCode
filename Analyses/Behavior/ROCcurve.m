function [Curves, Thresholds, AUC, StimIDs] = ROCcurve(ROIdata, ROIid, FrameIndex, TrialIndex, StimIndex, varargin)

type = 'dFoF'; % 'dFoF' or 'spikes'
Thresholds = -1:.01:15;

%% Parse input arguments
if ~exist('ROIdata','var') || isempty(ROIdata)
    directory = CanalSettings('DataDirectory');
    [ROIdata, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if ~exist('ROIid', 'var') || isempty(ROIid) || (ischar(ROIid) && strcmp(ROIid, 'all'))
    ROIid = [1 inf];
end

if ~exist('TrialIndex', 'var') || isempty(TrialIndex)
    TrialIndex = [1 inf];
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case {'Thresholds', 'thresholds', 'thresh'}
                Thresholds = varargin{index+1};
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

if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1),TrialIndex(1:end-1)+1:numel(ROIdata.DataInfo.TrialIndex)];
end
totalTrials = numel(TrialIndex);

if ~exist('StimIndex', 'var')
    StimIndex = ROIdata.DataInfo.StimID(TrialIndex);
end
StimIDs = unique(StimIndex);
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
numDataPoints = numel(Thresholds);

Curves = nan(numDataPoints, numStims, numROIs, 3);
AUC = nan(numROIs, numStims);
for rindex = 1:numROIs
    
    switch type
        case 'dFoF'
            data = ROIdata.rois(rindex).dFoF(TrialIndex, FrameIndex(1):FrameIndex(2));
        case 'spikes'
            data = ROIdata.rois(rindex).spike(TrialIndex, FrameIndex(1):FrameIndex(2));
    end
    
    for sindex = 1:numStims
        
%         [Curves(sindex,2,:,rindex),Curves(sindex,1,:,rindex),~,AUC(rindex,sindex)] = perfcurve(StimIndex, max(data,[],2), StimIndex(sindex), 'TVals', Thresholds);
        
        currentTrials = StimIndex == StimIDs(sindex);
        numTrials = sum(currentTrials);
        for tindex = 1:numDataPoints
            Curves(tindex, sindex, rindex, 1) = sum(any(data(currentTrials) >= Thresholds(tindex), 2))/numTrials;
            Curves(tindex, sindex, rindex, 2) = sum(any(data(~currentTrials) >= Thresholds(tindex), 2))/(totalTrials-numTrials);
            Curves(tindex, sindex, rindex, 3) = abs(cross([1,1,0]-[0,0,0],[Curves(tindex, sindex, rindex, 2), Curves(tindex, sindex, rindex, 1),0] - [0,0,0]))/abs([1,1,0] - [0,0,0]); % d-prime?
        end
        [~,order] = sort(Curves(:, sindex, rindex, 2));
        AUC(rindex, sindex) = trapz(Curves(order, sindex, rindex, 2), Curves(order, sindex, rindex, 1));
        
    end
    
%     figure; hold on;
%     plot(Curves(:, :, rindex, 2), Curves(:, :, rindex, 1));
%     legend(cellstr(num2str(StimIDs)), 'Location', 'SouthEast');
%     plot([0;1],[0;1], 'r--');
end