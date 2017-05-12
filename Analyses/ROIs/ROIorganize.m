function [ROIdata, series, ROIindex] = ROIorganize(ROIdata, AnalysisInfo, frames, ROIindex, varargin)

saveOut = false;
saveFile = '';

TrialIndex = [1 inf];
secondsBefore = 1.5;
secondsAfter = 4.5;
SeriesVariables = {}; % strings of fieldnames in 'frames' struct to extract to be trial-wise

directory = cd;

%% Parse input arguments
if ~exist('ROIdata','var') || isempty(ROIdata)
    [ROIdata, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if (~exist('AnalysisInfo','var') || isempty(AnalysisInfo)) && (~exist('frames', 'var') || isempty(frames))
    [ExperimentFile, p] = uigetfile({'*.mat'},'Choose Experiment file',directory);
    if isnumeric(ExperimentFile)
        return
    end
    ExperimentFile = fullfile(p,ExperimentFile);
end

if ~exist('ROIindex', 'var') || isempty(ROIindex)
    ROIindex = [1 inf];
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Trials','trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'secondsBefore'
                secondsBefore = varargin{index+1};
                index = index + 2;
            case 'secondsAfter'
                secondsAfter = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
                index = index + 2;
            case 'SeriesVariables'
                SeriesVariables = varargin{index+1};
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

fprintf('Organizing signals to be trial-wise...');

if ~isfield(ROIdata,'Config')
    warning('ROIdata does not have image ''Config'' struct attached; assuming numDepths=1 & frameRate=15.49');
    ROIdata.Config.Depth = 1;
    ROIdata.Config.FrameRate = 15.49;
    ROIdata.depth = 1;
end

%% Load stimulus information and determine trials to pull-out
if ischar(AnalysisInfo)
    ExperimentFile = AnalysisInfo;
    AnalysisInfo = [];
elseif ischar(frames)
    ExperimentFile = frames;
    frames = [];
end

if ~exist('AnalysisInfo', 'var') || isempty(AnalysisInfo)
    load(ExperimentFile, 'AnalysisInfo', '-mat');
end 

if ~exist('frames', 'var') || isempty(frames)
    load(ExperimentFile, 'frames', '-mat');
end

% Determine trials to extract
if TrialIndex(end)==inf
    TrialIndex = cat(2, TrialIndex(end-1), TrialIndex(end-1)+1:size(AnalysisInfo, 1));
end
TrialIndex(AnalysisInfo.ImgIndex(TrialIndex)==0) = []; % remove non-imaged trials
numTrials = numel(TrialIndex);

% Determine frame indices for each trial
numFramesBefore = round(secondsBefore*ROIdata.Config.FrameRate);
numFramesAfter = round(secondsAfter*ROIdata.Config.FrameRate);
depthID = idDepth(ROIdata.Config,[],'Depth',ROIdata.depth);
numFramesBefore = round(numFramesBefore/ROIdata.Config.Depth);
numFramesAfter = round(numFramesAfter/ROIdata.Config.Depth);

% OLD
% if isempty(numFramesBefore) % defaults to load all frames in the trial before the stimulus
%     numFramesBefore = mode(AnalysisInfo.TrialStimFrames(TrialIndex, 1)) - 1;
% end
% firstFrame = numFramesBefore + 1;
% if isempty(numFramesAfter) % defaults to load all frames in the trial after the stimulus starts
%     numFramesAfter = max(AnalysisInfo.nFrames(TrialIndex)) - firstFrame;
% end


%% Determine series data to reshape
if ~isempty(frames)
    seriesNames = fieldnames(frames);
    if ischar(SeriesVariables)
        SeriesVariables = {SeriesVariables};
    end
    for sindex = numel(SeriesVariables):-1:1
        if ~strcmp(SeriesVariables{sindex}, seriesNames)
            warning('Series variable ''%s'' not found, continuing without formatting it', SeriesVariables{sindex});
            SeriesVariables(sindex) = [];
        end
    end
else
    SeriesVariables = {};
end


%% Load ROI information and determine ROIs to reshape
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat'); % load in roi data
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end

% Determine ROIs to reshape
if ischar(ROIindex)
    switch ROIindex
        case {'all', 'All'}
            ROIindex = 1:numel(ROIdata.rois);
        case {'new', 'New'}
            ROIindex = find(arrayfun(@(x) (isempty(x.data)), ROIdata.rois));
    end
elseif isnumeric(ROIindex) && ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIindex);

% Organize neuropil data?
neuropil = false;
if isfield(ROIdata.rois, 'rawneuropil')
    neuropil = true;
end

% Organize spiking data?
spikes = false;
if isfield(ROIdata.rois, 'rawspikes')
    spikes = true;
end


%% Format each type of data to be trial-wise
ROIdata.DataInfo.TrialIndex = TrialIndex;
ROIdata.DataInfo.StimID = AnalysisInfo.StimID(TrialIndex);
ROIdata.DataInfo.numFramesBefore = numFramesBefore;
% ROIdata.DataInfo.numStimFrames = diff(AnalysisInfo.ExpStimFrames(TrialIndex,:),[],2) + 1; % OLD
ROIdata.DataInfo.numStimFrames = nan(numTrials,1);
ROIdata.DataInfo.numFramesAfter = numFramesAfter;

totalFrames = numFramesBefore + 1 + numFramesAfter;
FrameIndex = [AnalysisInfo.ExpStimFrames(TrialIndex, 1) - numFramesBefore, AnalysisInfo.ExpStimFrames(TrialIndex, 1) + numFramesAfter];

% Format series variables
if ~isempty(SeriesVariables)
    for sindex = 1:numel(SeriesVariables)
        series.(SeriesVariables{sindex}) = nan(numTrials, totalFrames);
        for nindex = 1:numTrials
            series.(SeriesVariables{sindex})(nindex,:) = frames.(SeriesVariables{sindex})(FrameIndex(nindex,1):FrameIndex(nindex,2));
        end
    end
else
    series = [];
end

% Format ROIs
for rindex = ROIindex
    ROIdata.rois(rindex).data = nan(numTrials, totalFrames);
    if neuropil
        ROIdata.rois(rindex).neuropil = nan(numTrials, totalFrames);
    end
    if spikes
        ROIdata.rois(rindex).spikes = nan(numTrials, totalFrames);
    end
    for nindex = 1:numTrials
        StimFrames = AnalysisInfo.ExpStimFrames(TrialIndex(nindex),:);
        relativeID = [find(depthID>=StimFrames(1),1,'first'), find(depthID<StimFrames(2),1,'last')]; % assigns stim frames to be frame during which stim came in through frame before frame during which stim moves out
        ROIdata.DataInfo.numStimFrames(nindex) = diff(relativeID) + 1;
        currentIndex = relativeID(1)-numFramesBefore:relativeID(1)+numFramesAfter;
        try
            ROIdata.rois(rindex).data(nindex,:) = ROIdata.rois(rindex).rawdata(currentIndex);
            if neuropil
                ROIdata.rois(rindex).neuropil(nindex,:) = ROIdata.rois(rindex).rawneuropil(currentIndex);
            end
            if spikes
                ROIdata.rois(rindex).spikes(nindex,:) = ROIdata.rois(rindex).rawspikes(currentIndex);
            end
            
        catch
            warning('More frames requested than in ROIdata, filling rest with NaNs');
            temp = ROIdata.rois(rindex).rawdata(currentIndex(1):end);
            ROIdata.rois(rindex).data(nindex,1:numel(temp)) = temp;
            if neuropil
                ROIdata.rois(rindex).neuropil(nindex,1:numel(temp)) = ROIdata.rois(rindex).rawneuropil(currentIndex(1):end);
            end
            if spikes
                ROIdata.rois(rindex).spikes(nindex,1:numel(temp)) = ROIdata.rois(rindex).rawspikes(currentIndex(1):end);
            end
        end
    end
end

fprintf('\tComplete\n');


%% Save to file
if saveOut
    % Save basic data
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', 'AnalysisInfo', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', 'AnalysisInfo', '-mat', '-append');
    end
    % Save series data
    if ~isempty(series)
        save(saveFile, 'series', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end
