function ProcessROIs(ImageFiles, ROIFiles, ExperimentFiles, varargin)

ExtractSignals = false;
EstimateSpikes = false;
OrganizeSignals = false;
ComputeDFoF = false;
ComputeAvgStim = false;

override = true; %index of files to force analysis of

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'ExtractSignals'}
                ExtractSignals = true;
                index = index + 1;
            case {'EstimateSpikes'}
                EstimateSpikes = true;
                index = index + 1;
            case {'OrganizeSignals'}
                OrganizeSignals = true;
                index = index + 1;
            case {'ComputeDFoF'}
                ComputeDFoF = true;
                index = index + 1;
            case {'ComputeAvgStim'}
                ComputeAvgStim = true;
                index = index + 1;
            case {'override'}
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
    directory = CanalSettings('DataDirectory');
    directory = cd;
    ImageFiles = uipickfiles('Prompt', 'Choose sbx files to process', 'FilterSpec', [directory, '*.sbx']);
    if isnumeric(ImageFiles)
        return
    end
end

if ~exist('ROIFiles', 'var') || isempty(ROIFiles)
    directory = CanalSettings('DataDirectory');
    ROIFiles = uipickfiles('Prompt', 'Choose ROI files to process', 'FilterSpec', [directory, '*.ciexp']);
    if isnumeric(ROIFiles)
        return
    end
end

if ~exist('ExperimentFiles', 'var') || isempty(ExperimentFiles)
    directory = CanalSettings('ExperimentDirectory');
    ExperimentFiles = uipickfiles('Prompt', 'Choose experiment files to process', 'FilterSpec', [directory, '*.ciexp']);
    if isnumeric(ExperimentFiles)
        return
    end
end

numFiles = numel(ImageFiles);

%% Determine what has already been accomplished in input Experiment Files
variables = cell(numFiles,1);
% for index = 1:numFiles
%     variables{index} = whos(matfile(ROIFiles{index}));
% end

%% Extract ROI signals
if ExtractSignals
    for index = 1:numFiles
        if (islogical(override) && override == true) || ismember(index,override) || ~any(strcmp({variables{index}.name}, 'ImageFile'))
            fprintf('\nFile %d of %d:\t', index, numFiles);
            extractSignals(ImageFiles{index}, ROIFiles{index}, 'all', 'Save', 'GPU', 'MotionCorrect', ExperimentFiles{index});
        end
    end
end

%% Estimate Spike Timing
if EstimateSpikes
    for index = 1:numFiles
        if (islogical(override) && override == true) || ismember(index,override) || ~any(strcmp({variables{index}.name}, 'Spikes'))
            fprintf('\nFile %d of %d:\t', index, numFiles);
            estimateSpikeTiming(ROIFiles{index});
        end
    end
end

%% Process ROI signals
if OrganizeSignals
    for index = 1:numFiles
        if (islogical(override) && override == true) || ismember(index,override) || ~any(strcmp({variables{index}.name}, 'AnalysisInfo'))
            fprintf('\nFile %d of %d:\t', index, numFiles);
            ROIorganize(ROIFiles{index}, ExperimentFiles{index});
        end
    end
end

%% Compute dF/F of signals
if ComputeDFoF
    for index = 1:numFiles
        fprintf('\nFile %d of %d:\t', index, numFiles);
        computeDFoF(ROIFiles{index});
    end
end

%% Fit ROI signals
if ComputeAvgStim
    for index = 1:numFiles
        fprintf('\nFile %d of %d:\t', index, numFiles);
        StimSignificance(ROIFiles{index},0);
    end
end