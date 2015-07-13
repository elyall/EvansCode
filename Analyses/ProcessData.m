function ProcessData(ImageFiles, ExperimentFiles, varargin)
warning('off','MATLAB:handle_graphics:exceptions:SceneNode'); %tex warning with waitbar

MotionCorrect = false;
ProcessFrames = false;
ProcessExperiment = false;
compAvgResponse = false;

override = false;

%% Process input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'MotionCorrect'}
                MotionCorrect = true;
                index = index + 1;
            case {'ProcessFrames'}
                ProcessFrames = true;
                index = index + 1;
            case {'ProcessExperiment'}
                ProcessExperiment = true;
                index = index + 1;
            case {'compAvgResponse'}
                compAvgResponse = true;
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
%     directory = CanalSettings('DataDirectory');
    directory = cd;
    ImageFiles = uipickfiles('Prompt', 'Choose sbx files to process', 'FilterSpec', [directory, '*.sbx']);
    if isnumeric(ImageFiles)
        return
    end
    [ImageFiles,p] = uigetfile({'*.sbx;*.tif;*.imgs'},'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ImageFiles)
        return
    elseif ischar(ImageFiles)
        ImageFiles = {fullfile(p, ImageFiles)};
    else
        for index = 1:numel(ImageFiles)
            ImageFiles{index} = fullfile(p,ImageFiles{index});
        end
    end
end

if ~exist('ExperimentFiles', 'var') || isempty(ExperimentFiles)
%     directory = CanalSettings('ExperimentDirectory');
    ExperimentFiles = uipickfiles('Prompt', 'Choose experiment files to process', 'FilterSpec', [directory, '*.ciexp']);
    if isnumeric(ExperimentFiles)
        return
    end
    [ExperimentFiles,p] = uigetfile({'*.ciexp;*.mat'},'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ExperimentFiles)
        return
    elseif ischar(ExperimentFiles)
        ExperimentFiles = {fullfile(p, ExperimentFiles)};
    else
        for index = 1:numel(ExperimentFiles)
            ExperimentFiles{index} = fullfile(p,ExperimentFiles{index});
        end
    end
end

%% Pull out all sbx files from any input directories
% FIX TO SELECT CORRESPONDING EXPERIMENT FILES
% nItems = length(ImageFiles);
% for I = nItems:-1:1
%     if isdir(ImageFiles{I})
%         files = dir([ImageFiles{I},'\*.sbx']);
%         files = {files(:).name};
%         for F = 1:length(files)
%             files{F} = fullfile(ImageFiles{I}, files{F});
%         end
%         ImageFiles(I) = [];
%         ImageFiles = [ImageFiles, files];
%     end
% end

%% Determine number of files
numFiles = numel(ImageFiles);

%% Determine what has already been accomplished in input Experiment Files
variables = cell(numFiles,1);
for index = 1:numFiles
    variables{index} = whos(matfile(ExperimentFiles{index}));
end

%% Post-process images
%Scanbox v1.0-1.2: Cycle through each file converting it to a .imgs file
% nFiles = length(SbxFiles);
% ImageFiles = cell(nFiles,1);
% for F = 1:nFiles
%     ImageFiles{F} = convertSbx(SbxFiles{F});
% end

%% Perform motion correction
if MotionCorrect
    for index = 1:numFiles;
        if ~any(strcmp({variables{index}.name}, 'MCdata')) || override
            fprintf('\nFile %d of %d:\t', index, numFiles);
            fullDoLucasKanade(ImageFiles{index}, [], 'SaveAlignmentTo', ExperimentFiles{index});
        end
    end
end

%% Compute projections
if ProcessFrames
    for index = 1:numFiles;
        if ~any(strcmp({variables{index}.name}, 'ImageFiles')) || override
            fprintf('\nFile %d of %d:\t', index, numFiles);
            computeProjections(ImageFiles{index}, [2,inf], ExperimentFiles{index}); %[2,500]
        end
    end
end

%% Determine Stimulus Frames
if ProcessExperiment
    for index = 1:numFiles;
        if ~any(strcmp({variables{index}.name}, 'AnalysisInfo')) || override
            fprintf('\nFile %d of %d:\t', index, numFiles);
            sbxPostProcess(ExperimentFiles{index}, ImageFiles{index}, true);
        end
    end
end

%% Compute Average Response for Each Stimulus & Create Tuning Curves
if compAvgResponse
    for index = 1:numFiles;
        if ~any(strcmp({variables{index}.name}, 'AvgEvokedDFoF')) || override
            fprintf('\nFile %d of %d:\t', index, numFiles);
            computeAverageStimResponse(ImageFiles{index}, ExperimentFiles{index});
        end
    end
end

%% Fit Tuning Curves for Vertical Bar Stimulus
% fitWholeFoV(ExperimentFile)

