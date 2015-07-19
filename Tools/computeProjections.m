function [AverageFrame, MaxProjection, MinProjection, Variance, Kurtosis] = computeProjections(Images, FrameIndex, ExperimentFile, varargin)

computeAverage = false;
computeMax = false;
computeMin = false;
computeVar = false; % dependent upon average
computeKur = false; % dependent upon average and variance

saveOut = false;
saveFile = ''; % defaults to ExperimentFile if input
loadType = 'Direct'; % 'MemMap' or 'Direct'
loadPrevious = false;
MotionCorrect = false; % false, filename to load MCdata from, or true to prompt for file selection

% default settings
framedim = 5;

% Memory settings
portionOfMemory = 0.08; % find 10% or less works best
sizeRAM = 32000000000; % amount of memory on your computer (UNIX-only)

%% Initialize Parameters
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Avg','avg'}
                if numel(varargin{index+1}) == 1
                    computeAverage = varargin{index+1};
                else
                    computeAverage = true;
                    AverageFrame = varargin{index+1};
                end
                index = index + 2;
            case {'Max','max'}
                if numel(varargin{index+1}) == 1
                    computeMax = varargin{index+1};
                else
                    computeMax = true;
                    MaxProjection = varargin{index+1};
                end
                index = index + 2;
            case {'Min','min'}
                if numel(varargin{index+1}) == 1
                    computeMin = varargin{index+1};
                else
                    computeMin = true;
                    MinProjection = varargin{index+1};
                end
                index = index + 2;
            case {'Var','var'}
                if numel(varargin{index+1}) == 1
                    computeVar = varargin{index+1};
                else
                    computeVar = true;
                    Variance = varargin{index+1};
                end
                index = index + 2;
            case {'Kur','kur'}
                if numel(varargin{index+1}) == 1
                    computeKur = varargin{index+1};
                else
                    computeKur = true;
                    Kurtosis = varargin{index+1};
                end
                index = index + 2;
            case 'MotionCorrect'
                MotionCorrect = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case 'SaveFile'
                saveFile = varargin{index+1};
                index = index + 2;
            case 'loadType'
                loadType = varargin{index+1};
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

if ~exist('Images', 'var') || isempty(Images)
    directory = CanalSettings('DataDirectory');
    [Images,p] = uigetfile({'*.imgs;*.sbx;*.tif'},'Select images:',directory,'MultiSelect','on');
    if isnumeric(Images)
        return
    elseif iscell(Images)
        for findex = 1:numel(Images)
            Images{findex} = fullfile(p, Images{findex});
        end
    elseif ischar(Images)
        Images = {fullfile(p, Images)};
    end
elseif ischar(Images)
    Images = {Images};
end

if ((saveOut && isempty(saveFile)) || isequal(MotionCorrect, true) || loadPrevious) && (~exist('ExperimentFile','var') || isempty(ExperimentFile))
    directory = CanalSettings('ExperimentDirectory');
    [ExperimentFile, p] = uigetfile({'*.mat'},'Choose Experiment file to load from and/or save to',directory);
    if isnumeric(ExperimentFile)
        if isequal(MotionCorrect, true)
            MotionCorrect = false;
        end
        if loadPrevious
            loadPrevious = false;
        end
        if saveOut && isempty(saveFile)
            saveOut = false;
        end
    else
        ExperimentFile = fullfile(p,ExperimentFile);
        if isequal(MotionCorrect, true)
            MotionCorrect = ExperimentFile;
        end
        if saveOut && isempty(saveFile)
            saveFile = ExperimentFile;
        end
    end
end

if ~exist('FrameIndex', 'var') || isempty(FrameIndex)
    FrameIndex = [1, inf];
end

%% Determine file to save to
if saveOut && isempty(saveFile)
    saveFile = ExperimentFile;
end

%% Load in Data and determine dimensions
if iscellstr(Images) % filename input
    ImageFile = Images;
    switch loadType
        case 'MemMap'
            [Images, loadObj, Config] = load2P(ImageFile, 'Type', 'MemMap', 'Double');
            if strcmp(loadObj.files(1).ext, '.sbx')
                Images2 = Images;
                clear Images;
                Images = intmax(loadObj.Precision) - Images2;
            end
            numFrames = sum([Config(:).Frames]);
            numFramesPerLoad = numFrames; % all frames already mapped
        case 'Direct'
            [Images, loadObj, Config] = load2P(ImageFile, 'Type', 'Direct', 'Frames', 2, 'Double');
            sizeFrame = whos('Images');
            sizeFrame = sizeFrame.bytes;
            if ispc
                mem = memory;
                numFramesPerLoad = max(1, floor(portionOfMemory*mem.MaxPossibleArrayBytes/sizeFrame));
            else
                numFramesPerLoad = max(1, floor(portionOfMemory*sizeRAM/sizeFrame));
            end
            numFrames = sum([Config(:).Frames]);
    end
    dim = loadObj.size;
    
else % numeric array input
    loadType = false;
    dim = size(Images);
    numFrames = size(Images, 5);
    numFramesPerLoad = numFrames; % all frames already loaded
    ImageFile = {''};
end
dim(framedim) = [];


%% Determine frames to analyze
if FrameIndex(end) == inf
    FrameIndex = cat(2, FrameIndex(1:end-1), FrameIndex(end-1)+1:numFrames);
end
totalFrames = numel(FrameIndex);


%% Load Previous Values
if loadPrevious
    load(ExperimentFile, 'ImageFiles');
    if exist('ImageFiles', 'var')
        if isfield(ImageFiles, 'Average')
            computeAverage = false;
            AverageFrame = ImageFiles.Average;
        end
        if isfield(ImageFiles, 'Min')
            computeMin = false;
        end
        if isfield(ImageFiles, 'Max');
            computeMax = false;
        end
        if isfield(ImageFiles, 'Var');
            computeVar = false;
            Variance = ImageFiles.Var;
        end
        if isfield(ImageFiles, 'Kur');
            computeKur = false;
        end
    end
end

%% Load in motion correction information
if ischar(MotionCorrect) % filename input
    load(MotionCorrect, 'MCdata', '-mat');
    MotionCorrect = true;
elseif isequal(MotionCorrect, true) % load from ExperimentFile
    load(ExperimentFile, 'MCdata', '-mat');
end
if ~exist('MCdata', 'var')
    MotionCorrect = false;
end

%% Compute Average
if computeAverage || computeMax || computeMin
    
    % Initialize outputs
    if computeAverage && ~exist('AverageFrame', 'var')
        AverageFrame = zeros(dim);
    end
    if computeMax && ~exist('MaxProjection', 'var')
        MaxProjection = zeros(dim);
    end
    if computeMin && ~exist('MinProjection', 'var')
        MinProjection = inf(dim);
    end
    
    % Cycle through frames computing projections
    fprintf('Computing projections of %d frames: %s\n', totalFrames, ImageFile{1}); % update status
    for bindex = 1:numFramesPerLoad:totalFrames % direct loading only -> load frames in batches
        lastframe = min(bindex+numFramesPerLoad-1, totalFrames);
        currentFrames = FrameIndex(bindex:lastframe);
        
        % direct loading only -> load current batch
        if strcmp(loadType, 'Direct')
            [Images, loadObj] = load2P(ImageFile, 'Type', 'Direct', 'Frames', currentFrames, 'Double'); %direct
        end
        
        % Correct for motion
        if MotionCorrect
            Images = applyMotionCorrection(Images, MCdata, loadObj);
        end
        
        % Compute average frame
        if computeAverage
            AverageFrame = AverageFrame + sum(Images, framedim)/totalFrames;
        end
        
        % Compute maximum projection
        if computeMax
            MaxProjection = max(cat(framedim, MaxProjection, Images), [], framedim);
        end
        
        % Compute minimum projection
        if computeMin
            MinProjection = min(cat(framedim, MinProjection, Images), [], framedim);
        end
        
        % Compute variance and kurtosis
        if numFramesPerLoad >= FrameIndex(2)
            if computeVar
                meanSubtracted = bsxfun(@minus, Images, AverageFrame);
                Variance = mean(meanSubtracted.^2, framedim); % should be equivalent to 'var(Images,0,5)'
            end
            if computeKur
                Kurtosis = mean(bsxfun(@rdivide, meanSubtracted, sqrt(Variance)).^4, framedim) - 3; % should be equivalent to 'kurtosis(Images,1,5)'
            end
        end
        
        fprintf('\tFinished frames %d through %d\n', bindex, lastframe); % display update
    end
end

% DIRECT ONLY: Can't load all images at once => compute variance
% Variance calculation requires knowing the mean
if numFramesPerLoad < FrameIndex(2) && computeVar
    fprintf('Computing Variance...\n');
    
    if ~exist('Variance', 'var')
        Variance = zeros(dim);
    end
    
    for bindex = 1:numFramesPerLoad:totalFrames % direct loading only -> load frames in batches
        lastframe = min(bindex+numFramesPerLoad-1, totalFrames);
        currentFrames = FrameIndex(bindex:lastframe);
        
        % direct loading only -> load current batch
        if strcmp(loadType, 'Direct')
            [Images, loadObj] = load2P(ImageFile, 'Type', 'Direct', 'Frames', currentFrames, 'Double'); %direct
        end
        
        % Correct for motion
        if MotionCorrect
            Images = applyMotionCorrection(Images, MCdata, loadObj);
        end
        
        % Compute variance
        Variance = Variance + sum(bsxfun(@minus, Images, AverageFrame).^2, framedim)/totalFrames; %var = 1/T*sum((x_t-x_avg)^2)
        
        fprintf('\tVariance: Finished frames %d through %d\n', bindex, lastframe); % display update
    end
end

% DIRECT ONLY: Can't load all images at once => compute kurtosis
% Kurtosis calculation requires knowing the mean and variance
if numFramesPerLoad < FrameIndex(2) && computeKur
    fprintf('Computing Kurtosis...\n');
    
    if ~exist('Kurtosis', 'var')
        Kurtosis = zeros(dim);
    end
    
    for bindex = 1:numFramesPerLoad:totalFrames % direct loading only -> load frames in batches
        lastframe = min(bindex+numFramesPerLoad-1, totalFrames);
        currentFrames = FrameIndex(bindex:lastframe);
        
        % direct loading only -> load current batch
        if strcmp(loadType, 'Direct')
            [Images, loadObj] = load2P(ImageFile, 'Type', 'Direct', 'Frames', currentFrames, 'Double'); %direct
        end
        
        % Correct for motion
        if MotionCorrect
            Images = applyMotionCorrection(Images, MCdata, loadObj);
        end
        
        % Compute kurtosis
        Kurtosis = Kurtosis + sum(bsxfun(@rdivide, bsxfun(@minus, Images, AverageFrame), sqrt(Variance)).^4, framedim)/fileConfig.Frames - 3; %kur = 1/T*sum(((x_t-x_avg)/x_var)^4)
        
        fprintf('\tKurtosis: Finished frames %d through %d\n', bindex, lastframe); % display update
    end
end

%% Save to file
if saveOut
    % Append to file
    %     load(ExperimentFile, 'ImageFiles', '-mat');
    %     if ~exist('ImageFiles', 'var')
    index = 1;
    %     else
    %         if any(strcmp({ImageFiles(:).filename}, ImgsFile));
    %             index = find(strcmp({ImageFiles(:).filename}, ImgsFile)); % replace
    %         else
    %             index = numel(ImageFiles + 1); % add to end
    %         end
    %     end
    ImageFiles(index).filename = ImageFile;
    ImageFiles(index).frames = FrameIndex;
    if computeAverage
        ImageFiles(index).Average = AverageFrame;
    end
    if computeMax
        ImageFiles(index).Max = MaxProjection;
    end
    if computeMin
        ImageFiles(index).Min = MinProjection;
    end
    if computeVar
        ImageFiles(index).Var = Variance;
    end
    if computeKur
        ImageFiles(index).Kur = Kurtosis;
    end
    if ~exist(saveFile, 'file')
        save(saveFile, 'ImageFiles', '-mat');
    else
        save(saveFile, 'ImageFiles', '-mat', '-append');
    end
    fprintf('Saved projections to: %s\n', saveFile);
end