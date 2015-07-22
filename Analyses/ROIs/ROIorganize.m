function [ROIdata, series, ROIid] = ROIorganize(ROIdata, AnalysisInfo, frames, ROIid, varargin)

saveOut = false;

TrialIndex = [1 inf];
numFramesBefore = []; % '[]' for default
numFramesAfter = []; % '[]' for default
saveFile = '';
SeriesVariables = {}; % strings of fieldnames in 'frames' struct to extract to be trial-wise


%% Parse input arguments
if ~exist('ROIdata','var') || isempty(ROIdata)
    directory = CanalSettings('DataDirectory');
    [ROIdata, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if (~exist('AnalysisInfo','var') || isempty(AnalysisInfo)) && (~exist('frames', 'var') || isempty(frames))
    directory = CanalSettings('ExperimentDirectory');
    [ExperimentFile, p] = uigetfile({'*.mat'},'Choose Experiment file',directory);
    if isnumeric(ExperimentFile)
        return
    end
    ExperimentFile = fullfile(p,ExperimentFile);
end

if ~exist('ROIid', 'var') || isempty(ROIid)
    ROIid = 'all'; % 'all' or vector of indices
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Trials','trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'numFramesBefore'
                numFramesBefore = varargin{index+1};
                index = index + 2;
            case 'numFramesAfter'
                numFramesAfter = varargin{index+1};
                index = index + 2;
            case 'Save'
                saveOut = true;
                index = index + 1;
            case 'SaveFile'
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
if isempty(numFramesBefore) % defaults to load all frames in the trial before the stimulus
    numFramesBefore = mode(AnalysisInfo.TrialStimFrames(TrialIndex, 1)) - 1;
end
firstFrame = numFramesBefore + 1;
if isempty(numFramesAfter) % defaults to load all frames in the trial after the stimulus starts
    numFramesAfter = max(AnalysisInfo.nFrames(TrialIndex)) - firstFrame;
end
totalFrames = numFramesBefore + 1 + numFramesAfter;
FrameIndex = [AnalysisInfo.ExpStimFrames(TrialIndex, 1) - numFramesBefore, AnalysisInfo.ExpStimFrames(TrialIndex, 1) + numFramesAfter];


%% Determine series data to reshape
seriesNames = fieldnames(frames);
for sindex = numel(SeriesVariables):-1:1
    if ~strcmp(SeriesVariables{sindex}, seriesNames)
        warning('Series variable %s not found', SeriesVariables{sindex});
        SeriesVariables(sindex) = [];
    end
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
numROIs = numel(ROIdata.rois);

% Determine ROIs to extract signals for
if ischar(ROIid) && strcmp(ROIid, 'all')
    numROIs = numel(ROIdata.rois);
    ROIid = 1:numROIs;
end

% Organize neuropil data?
neuropil = false;
if isfield(ROIdata.rois, 'rawneuropil')
    neuropil = true;
end

% Organize spiking data
spikes = false;
if isfield(ROIdata.rois, 'rawspikes') % not sure if this is correct, check upstream function that estimates spikes to see output
    spikes = true;
end

% Fix fields (for older versions)
if ~isfield(ROIdata.rois, 'rawdata')
    for rindex = 1:numROIs
        ROIdata.rois(rindex).rawdata = ROIdata.rois(rindex).data;
        if neuropil
            ROIdata.rois(rindex).rawneuropil = ROIdata.rois(rindex).neuropil;
        end
    end
end


%% Extract ROI data to one matrix
numFrames = numel(ROIdata.rois(1).rawdata);
Data = zeros(numROIs, numFrames);
if neuropil
    Neuropil = zeros(numROIs, numFrames);
end
% if spikes
%     Spikes = zeros(numROIs, numFrames);
% end
for rindex = 1:numROIs
    Data(rindex,:) = ROIdata.rois(rindex).rawdata;
    if neuropil
        Neuropil(rindex,:) = ROIdata.rois(rindex).rawneuropil;
    end
%     if spikes
%         Spikes(rindex,:) = ROIdata.rois(rindex).rawspikes;
%     end
end


%% Format each type of data to be trial-wise
ROIdata.DataInfo.TrialIndex = TrialIndex;
ROIdata.DataInfo.StimID = AnalysisInfo.StimID(TrialIndex);
ROIdata.DataInfo.numFramesBefore = numFramesBefore;
ROIdata.DataInfo.numStimFrames = diff(AnalysisInfo.TrialStimFrames(TrialIndex,:),[],2) + 1;
ROIdata.DataInfo.numFramesAfter = numFramesAfter;

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
for rindex = ROIid
    ROIdata.rois(rindex).data = nan(numTrials, totalFrames);
    if neuropil
        ROIdata.rois(rindex).neuropil = nan(numTrials, totalFrames);
    end
    if spikes
        ROIdata.rois(rindex).spikes = nan(numTrials, totalFrames);
    end
    for nindex = 1:numTrials
        try
            ROIdata.rois(rindex).data(nindex,:) = ROIdata.rois(rindex).rawdata(FrameIndex(nindex,1):FrameIndex(nindex,2));
            if neuropil
                ROIdata.rois(rindex).neuropil(nindex,:) = ROIdata.rois(rindex).rawneuropil(FrameIndex(nindex,1):FrameIndex(nindex,2));
            end
            if spikes
                ROIdata.rois(rindex).spikes(nindex,:) = ROIdata.rois(rindex).rawspikes(FrameIndex(nindex,1):FrameIndex(nindex,2));
            end
            
        catch
            warning('More frames requested than in ROIdata, filling rest with NaNs');
            temp = ROIdata.rois(rindex).rawdata(FrameIndex(nindex,1):end);
            ROIdata.rois(rindex).data(nindex,1:numel(temp)) = temp;
            if neuropil
                ROIdata.rois(rindex).neuropil(nindex,1:numel(temp)) = ROIdata.rois(rindex).rawneuropil(FrameIndex(nindex,1):end);
            end
            if spikes
                ROIdata.rois(rindex).spikes(nindex,1:numel(temp)) = ROIdata.rois(rindex).rawspikes(FrameIndex(nindex,1):end);
            end
        end
    end
end


%% Save to file
if saveOut
    % Save basic data
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', 'Data', 'Neuropil', 'AnalysisInfo');
    else
        save(saveFile, 'ROIdata', 'Data', 'Neuropil', 'AnalysisInfo', '-append');
    end
    % Save series data
    if ~isempty(series)
        save(saveFile, 'series', '-append');
    end
    fprintf('ROIdata saved to: %s\n', saveFile);
end
