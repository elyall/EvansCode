function out = computeTrialMean2(Data,FrameIndex,Fieldname)
% Data is  numFrames x numTrials x numROIs
% out is numROIs x numTrials

directory = cd;


%% Check input arguments
if ~exist('Data','var') || isempty(Data)
    [Data, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file',directory);
    if isnumeric(Data)
        return
    end
    Data = fullfile(p,Data);
end

if ~exist('FrameIndex','var') || isempty(FrameIndex)
    FrameIndex = 6:10; % index of frames containing stimulus
end

if ~exist('Fieldname','var') || isempty(Fieldname)
    Fieldname = 'dFoF'; % 'dFoF' or 'spikes'
end


%% Load ROI data
if ischar(Data)
    load(Data, 'ROIdata', '-mat');
    Data = ROIdata;
end
if isstruct(Data)
    Data = gatherROIdata(Data,Fieldname);
end
numTrials = size(Data,2);
if isnumeric(FrameIndex)
    FrameIndex = {FrameIndex};
end
if numel(FrameIndex)==1 && numTrials~=1
    FrameIndex = repmat(FrameIndex,numTrials,1);
end


%% Calculate mean per trial

% out = squeeze(nanmean(Data(FrameIndex,:,:),1)); % nan frames may exist due to motion artifacts
out = arrayfun(@(x) squeeze(nanmean(Data(FrameIndex{x},x,:),1)), 1:numTrials,'UniformOutput',false);
out = cat(2,out{:});

