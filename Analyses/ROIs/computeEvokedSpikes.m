function [Data,Baseline] = computeEvokedSpikes(Data, FrameIndex)
% Data is  numFrames x numTrials x numROIs

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
    FrameIndex = 1:5; % index of baseline period
end


%% Load ROI data
if ischar(Data)
    load(Data, 'ROIdata', '-mat');
    Data = ROIdata;
end
if isstruct(Data)
    Data = gatherROIdata(Data,'rawspikes');
end


%% Compute evoked spikes for each trial for each ROI
fprintf('Calculating evoked spikes...');

Baseline = nanmean(Data(FrameIndex, :, :), 1); % avg across time for each trial
Baseline = nanmean(Baseline, 2); % avg across trials
Data = Data - Baseline;         % subtract off baseline spikes
Baseline = squeeze(Baseline);

fprintf('\tComplete\n');

