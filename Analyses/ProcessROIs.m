function ProcessROIs(ImageFiles, ROIFiles, ExperimentFiles, varargin)

% Average over ROIs
ExtractSignals = false;

% Event detection
EstimateSpikes = false;

% Trial-wise sorting
OrganizeSignals = false;

% Trial-wise dF/F
ComputeDFoF = false;

% Tuning Curves
ComputeAvgStim = false;
minrunspeed = 100;
outlierweight = 4;

directory = cd;
override = false; 

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ExtractSignals'
                ExtractSignals = true;
                index = index + 1;
            case 'EstimateSpikes'
                EstimateSpikes = true;
                index = index + 1;
            case 'OrganizeSignals'
                OrganizeSignals = true;
                index = index + 1;
            case 'ComputeDFoF'
                ComputeDFoF = true;
                index = index + 1;
            case 'ComputeAvgStim'
                ComputeAvgStim = true;
                index = index + 1;
            case 'minrunspeed'
                minrunspeed = varargin{index+1};
                index = index + 2;
            case 'outlierweight'
                outlierweight = varargin{index+1};
                index = index + 2;
            case 'override'
                override = true;
                index = index + 1;
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

%% Post-process ROIs

numFiles = numel(ImageFiles);
for findex = 1:numFiles;
    fprintf('Processing ROI files %d of %d: %s\n', findex, numFiles, ROIFiles{findex});
    
    
    %% Determine file to save to
    saveFile = ROIFiles{findex};
    
    
    %% Determine what has already been accomplished
    variables = whos(matfile(saveFile));
    
    
    %% Extract ROI signals
    if ExtractSignals && (~any(strcmp({variables.name}, 'ImageFiles')) || override)
        ROIdata = extractSignals(ImageFiles{findex}, ROIFiles{findex}, 'all', 'Save', 'SaveFile', saveFile, 'GPU', 'MotionCorrect', ExperimentFiles{findex});
    end
    
    
    %% Estimate Spike Timing
    if EstimateSpikes && (~any(strcmp({variables.name}, 'Spikes')) || override)
        if ~exist('ROIdata', 'var')
            ROIdata = estimateSpikeTiming(ROIFiles{findex}, 0.65, 'Save', 'SaveFile', saveFile);
        else
            ROIdata = estimateSpikeTiming(ROIFiles{findex}, 0.65, 'Save', 'SaveFile', saveFile);
        end
    end
    
    
    %% Sort ROI signals to be trial-wise
    if OrganizeSignals && (~any(strcmp({variables.name}, 'AnalysisInfo')) || override)
        if ~exist('ROIdata', 'var')
            ROIdata = ROIorganize(ROIFiles{findex}, ExperimentFiles{findex}, [], 'all', 'Save', 'SaveFile', saveFile);
        else
            ROIdata = ROIorganize(ROIFiles{findex}, ExperimentFiles{findex}, [], 'all', 'Save', 'SaveFile', saveFile);
        end
    end
    
    
    %% Load ROIdata
    if ~exist('ROIdata', 'var')
        load(ROIFiles{findex}, 'ROIdata', '-mat');
    end
    
    
    %% Compute dF/F of trials
    if ComputeDFoF && (~isfield(ROIdata.rois, 'dFoF') || override)
        ROIdata = computeDFoF(ROIdata, 0.65, 'Save', 'SaveFile', saveFile);
    end
    
    
    %% Compute tuning curves
    if ComputeAvgStim && (~isfield(ROIdata.rois, 'curve') || override)
        
        % Determine trials to analyze
        load(ExperimentFiles{findex}, 'AnalysisInfo', 'frames', '-mat');
        if exist('frames', 'var')
            TrialIndex = determineRunning(AnalysisInfo, frames, minrunspeed);
        else
            TrialIndex = [1 inf];
        end
        
        % Compute tuning curves
        computeTuningCurve(ROIdata, [1 inf], TrialIndex,...
            'Fit',...
            'ControlID', 0,...
            'outlierweight', outlierweight,...
            'Save',...
            'SaveFile', saveFile);
    end
    
    clear ROIdata
    
end %files