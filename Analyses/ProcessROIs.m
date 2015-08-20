function saveFiles = ProcessROIs(ImageFiles, ROIFiles, ExperimentFiles, varargin)

% Average over ROIs
ExtractSignals = false;
FrameIndex = {[1 inf]};

% Event detection
EstimateSpikes = false;

% Trial-wise sorting
OrganizeSignals = false;

% Trial-wise dF/F
ComputeDFoF = false;
NeuropilWeight = {[]};

% Tuning Curves
ComputeAvgStim = false;
minrunspeed = 100;
outlierweight = 4;

saveOut = true;
saveFiles = {[]};
ROIindex = [1 inf];
directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ExtractSignals'
                ExtractSignals = true;
                index = index + 1;
            case {'Frames', 'frames', 'FrameIndex'}
                FrameIndex = varargin{index+1};
                index = index + 2;
            case 'EstimateSpikes'
                EstimateSpikes = true;
                index = index + 1;
            case 'OrganizeSignals'
                OrganizeSignals = true;
                index = index + 1;
            case 'ComputeDFoF'
                ComputeDFoF = true;
                index = index + 1;
            case 'NeuropilWeight'
                NeuropilWeight = varargin{index+1};
                index = index + 2;
            case 'ComputeAvgStim'
                ComputeAvgStim = true;
                index = index + 1;
            case 'minrunspeed'
                minrunspeed = varargin{index+1};
                index = index + 2;
            case 'outlierweight'
                outlierweight = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFiles', 'saveFiles'}
                saveFiles = varargin{index+1};
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

if ~exist('ImageFiles', 'var') || isempty(ImageFiles)
    [ImageFiles,p] = uigetfile({'*.sbx;*.tif;*.imgs'}, 'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ImageFiles) % no file selected
        return
    elseif iscellstr(ImageFiles) % multiple files selected
        ImageFiles = fullfile(p, ImageFiles);
    elseif ischar(ImageFiles) % single file selected
        ImageFiles = {fullfile(p, ImageFiles)};
    end
elseif ischar(ImageFiles)
    if isdir(ImageFiles) % directory input
        p = ImageFiles;
        ImageFiles = dir(p);
        ImageFiles = fullfile(p, {ImageFiles(~cellfun(@isempty, regexpi({ImageFiles.name}, '.*(sbx|tif)'))).name});
    else % single file input
        ImageFiles = {ImageFiles};
    end
end

if ~exist('ROIFiles', 'var') || isempty(ROIFiles)
    [ROIFiles,p] = uigetfile({'*.rois;*.mat'}, 'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ROIFiles) % no file selected
        return
    elseif iscellstr(ROIFiles) % multiple files selected
        ROIFiles = fullfile(p, ROIFiles);
    elseif ischar(ROIFiles) % single file selected
        ROIFiles = {fullfile(p, ROIFiles)};
    end
elseif ischar(ROIFiles)
    if isdir(ROIFiles) % directory input
        p = ROIFiles;
        ROIFiles = dir(fullfile(p, '*.rois'));
        ROIFiles = fullfile(p, {ROIFiles.name});
    else % single file input
        ROIFiles = {ROIFiles};
    end
end

if ~exist('ExperimentFiles', 'var') || isempty(ExperimentFiles) % nothing input
    [ExperimentFiles,p] = uigetfile({'*.exp;*.mat'},'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ExperimentFiles) % no file selected
        return
    elseif iscellstr(ExperimentFiles) % multiple files selected
        ExperimentFiles = fullfile(p, ExperimentFiles);
    elseif ischar(ExperimentFiles) % single file selected
        ExperimentFiles = {fullfile(p, ExperimentFiles)};
    end
elseif ischar(ExperimentFiles)
    if isdir(ExperimentFiles) % directory input
        p = ExperimentFiles;
        ExperimentFiles = dir(fullfile(p, '*.exp'));
        ExperimentFiles = fullfile(p, {ExperimentFiles.name});
    else % single file input
        ExperimentFiles = {ExperimentFiles};
    end
end


%% Prepare indices & other independent variables
numFiles = numel(ROIFiles);

if ~iscell(ROIindex)
    ROIindex = {ROIindex};
end
if numel(ROIindex) == 1
    ROIindex = repmat(ROIindex, numFiles, 1);
end

if ~iscell(FrameIndex)
    FrameIndex = {FrameIndex};
end
if numel(FrameIndex) == 1
    FrameIndex = repmat(FrameIndex, numFiles, 1);
end

if ~iscell(NeuropilWeight)
    NeuropilWeight = {NeuropilWeight};
end
if numel(NeuropilWeight) == 1
    NeuropilWeight = repmat(NeuropilWeight, numFiles, 1);
end

if ~iscell(saveFiles)
    saveFiles = {saveFiles};
end
if numel(saveFiles) == 1
    saveFiles = repmat(saveFiles, numFiles, 1);
end


%% Cycle through files
for findex = 1:numFiles;
    fprintf('Processing ROI files %d of %d: %s\n', findex, numFiles, ROIFiles{findex});
    
    %% Determine file to save to
    if isempty(saveFiles{findex})
        saveFiles{findex} = ROIFiles{findex};
    end
    
    
    %% Determine what has already been accomplished
    variables = whos(matfile(saveFiles{findex}));
    
    
    %% Load ROIdata
    load(ROIFiles{findex}, 'ROIdata', '-mat');
    
    
    %% Extract ROI signals
    if ExtractSignals
        [ROIdata, Data, Neuropil, ROIindex{findex}] = extractSignals(ImageFiles{findex}, ROIdata, ROIindex{findex}, 'FrameIndex', FrameIndex{findex}, 'MotionCorrect', ExperimentFiles{findex});
        if saveOut && isequal(ROIindex{findex}, 1:numel(ROIdata.rois))
            save(saveFiles{findex}, 'Data', 'Neuropil', '-mat', '-append');
        end
    end
    
    if iscolumn(ROIindex{findex})
        ROIindex{findex} = ROIindex{findex}';
    end
    
    %% Estimate Spike Timing
    if EstimateSpikes
        [ROIdata, Spikes] = estimateSpikeTiming(ROIFiles{findex}, NeuropilWeight{findex}, 'ROIindex', ROIindex{findex});
        if saveOut && isequal(ROIindex{findex}, 1:numel(ROIdata.rois))
            save(saveFiles{findex}, 'Spikes', '-mat', '-append');
        end
    end
    
    
    %% Sort ROI signals to be trial-wise
    if OrganizeSignals
        [ROIdata, series] = ROIorganize(ROIdata, ExperimentFiles{findex}, [], ROIindex{findex}, 'SeriesVariables', 'RunningSpeed');
        if saveOut
            save(saveFiles{findex}, 'series', '-mat', '-append');
        end
    end
    
    
    %% Compute dF/F of trials
    if ComputeDFoF
        ROIdata = computeDFoF(ROIdata, 'ROIindex', ROIindex{findex}, 'NeuropilWeight', NeuropilWeight{findex});
    end
    
    
    %% Compute tuning curves
    if ComputeAvgStim
        
        % Determine trials to analyze
        load(ExperimentFiles{findex}, 'AnalysisInfo', 'frames', '-mat');
        if exist('frames', 'var') && isfield(frames, 'RunningSpeed')
            TrialIndex = determineRunning(AnalysisInfo, frames, minrunspeed);
        else
            TrialIndex = [1 inf];
        end
        
        % Compute tuning curves
        ROIdata = computeTuningCurve(ROIdata, ROIindex{findex}, TrialIndex,...
            'Fit',...
            'ControlID', 0,...
            'outlierweight', outlierweight);
    end
    
    
    %% Save data to file
    if saveOut
        if ~exist(saveFiles{findex}, 'file')
            save(saveFiles{findex}, 'ROIdata', '-mat', '-v7.3');
        else
            save(saveFiles{findex}, 'ROIdata', '-mat', '-append');
        end
        fprintf('\tROIdata saved to: %s\n', saveFiles{findex});
    end
    
    
end %files