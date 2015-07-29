function ProcessData(ImageFiles, ExperimentFiles, varargin)
warning('off','MATLAB:handle_graphics:exceptions:SceneNode'); %tex warning with waitbar

% Experiment
ProcessExperiment = false;

% Motion Correction
MotionCorrect = false;
numFramesTemplate = 500;

% Projections
ProcessFrames = false;

% Average Trial
compAvgResponse = false;
minrunspeed = 100;

directory = cd;
override = false;

%% Process input arguments
findex = 1;
while findex<=length(varargin)
    try
        switch varargin{findex}
            case 'MotionCorrect'
                MotionCorrect = true;
                findex = findex + 1;
            case 'ProcessFrames'
                ProcessFrames = true;
                findex = findex + 1;
            case 'ProcessExperiment'
                ProcessExperiment = true;
                findex = findex + 1;
            case 'compAvgResponse'
                compAvgResponse = true;
                findex = findex + 1;
            case 'override'
                override = true;
                findex = findex + 1;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{findex});
                findex = findex + 1;
        end
    catch
        warning('Argument %d not recognized',findex);
        findex = findex + 1;
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


%% Post-process images
%Scanbox v1.0-1.2: Cycle through each file converting it to a .imgs file
% nFiles = length(SbxFiles);
% ImageFiles = cell(nFiles,1);
% for F = 1:nFiles
%     ImageFiles{F} = convertSbx(SbxFiles{F});
% end

numFiles = numel(ImageFiles);
for findex = 1:numFiles;
    fprintf('Processing experiment %d of %d\n', findex, numFiles);
    RunningData = false;
    
    %% Determine file to save to
    if iscellstr(ExperimentFiles{findex})
        saveFile = ExperimentFiles{findex}{1};
    elseif ischar(ExperimentFiles{findex})
        saveFile = ExperimentFiles{findex};
    end
    
    
    %% Determine what has already been accomplished
    variables = whos(matfile(saveFile));

    
    %% Determine Stimulus Frames
    if ProcessExperiment && (~any(strcmp({variables.name}, 'AnalysisInfo')) || override)
        
        % Process experiment
        [~,frames] = sbxPostProcess2(ExperimentFiles{findex}, ImageFiles{findex}, 'Save', 'SaveFile', saveFile);
        
        % Determine if mouse's running speed was recorded
        if isfield(frames, 'RunningSpeed')
            RunningData = true;
        end
    elseif any(strcmp({variables.name}, 'frames')) % load running data for subsequent analyses
        load(ExperimentFiles{findex}, 'frames', '-mat');
        if isfield(frames, 'RunningSpeed')
            RunningData = true;
        end
    end
    
    
    %% Perform motion correction
    if MotionCorrect && (~any(strcmp({variables.name}, 'MCdata')) || override)
        
        % Make template
        if RunningData
            [~,TemplateIndices] = sort(frames.RunningSpeed);
            TemplateIndices(TemplateIndices==1) = []; % remove first frame (frame is incomplete in scanbox files)
            TemplateIndices = TemplateIndices(1:numFramesTemplate);
        else
            TemplateIndices = 2:2+numFramesTemplate-1;
        end
        
        % Perform motion correction
        fullDoLucasKanade(ImageFiles{findex}, TemplateIndices, 'SaveAlignmentTo', saveFile);
        
    end
    
    
    %% Compute projections
    if ProcessFrames && (~any(strcmp({variables.name}, 'ImageFiles')) || override)
        computeProjections(ImageFiles{findex}, [2,inf], saveFile, 'MotionCorrect', true, 'Avg', true, 'Min', true, 'Max', true, 'Var', true, 'Save', 'SaveFile', saveFile);
    end
    
    
    %% Compute Average Response for Each Stimulus & Create Tuning Curves
    if compAvgResponse && (~any(strcmp({variables.name}, 'AvgTrial')) || override)
        
        % Determine trials to analyze
        if RunningData
            TrialIndex = determineRunning(saveFile, [], minrunspeed);
        else
            TrialIndex = [1 inf];
        end
        
        % Compute average trial for each stimulus
        computeAverageStimResponse(ImageFiles{findex}, saveFile, TrialIndex, true, 'Save');
        
    end
    
end


