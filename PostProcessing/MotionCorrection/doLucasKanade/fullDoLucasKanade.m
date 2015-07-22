function [MCdata, Template] = fullDoLucasKanade(Images, Template, varargin)

% Default settings
Parameters.Nbasis = 16;
Parameters.niter = 30;
Parameters.damping = 1;
Parameters.deltacorr = .0009;

Channel2AlignFrom = 1;
numFramesInitialTemplate = 500;
numFramesSecondPass = 0;
MCdataFilename = false;
MCImgsFilename = false;
FrameIndex = [1 inf];

% Memory settings
portionOfMemory = 0.08; % find 10% or less works best
sizeRAM = 32000000000; % amount of memory on your computer (UNIX-only)

%% Parse input arguments
index = 1;
while index<=length(varargin)
    switch varargin{index}
        case 'SaveAlignmentTo'
            MCdataFilename = varargin{index+1};
            index = index + 2;
        case 'SaveImagesTo'
            MCImgsFilename = varargin{index+1};
            index = index + 2;
        case 'Nbasis'
            Parameters.Nbasis = varargin{index+1};
            index = index + 2;
        case 'niter'
            Parameters.niter = varargin{index+1};
            index = index + 1;
        case 'RefChannel'
            Channel2AlignFrom = varargin{index+1};
            index = index + 2;
        case 'MCdata'
            MCdata = varargin{index+1};
            index = index + 2;
        case 'Frames'
            FrameIndex = varargin{index+1};
            index = index + 2;
        otherwise
            warning('Argument ''%s'' not recognized',varargin{index});
            index = index + 1;
    end
end

if ~exist('Images','var') || isempty(Images) % Prompt for file selection
    directory = CanalSettings('DataDirectory');
    [ImageFiles,p] = uigetfile({'*.imgs;*.sbx;*.tif'},'Select images:',directory,'MultiSelect','on');
    if isnumeric(ImageFiles)
        return
    elseif iscell(ImageFiles)
        for findex = 1:numel(ImageFiles)
            ImageFiles{findex} = fullfile(p, ImageFiles{findex});
        end
    elseif ischar(ImageFiles)
        ImageFiles = {fullfile(p, ImageFiles)};
    end
    LoadImgs = true;
elseif ischar(Images) % File is specified
    ImageFiles = {Images};
    LoadImgs = true;
elseif iscellstr(Images)
    ImageFiles = Images;
    LoadImgs = true;
elseif ~isa(Images, 'double')
    Images = double(Images);
end

if ~exist('Template', 'var') || isempty(Template)
    Template = false;
end

%% Determine filenames to save images to
if MCdataFilename == true
    directory = CanalSettings('ExperimentDirectory');
    [MCdataFilename,p] = uigetfile({'*.mat'}, 'Save alignment data to:', directory);
    if isnumeric(MCdataFilename)
        return
    end
    MCdataFilename = fullfile(p, MCdataFilename);
end
if MCImgsFilename == true
    directory = CanalSettings('DataDirectory');
    [MCImgsFilename,p] = uiputfile({'*.sbx'}, 'Save registered images as:', directory);
    if isnumeric(MCImgsFilename)
        return
    end
    MCImgsFilename = fullfile(p, MCImgsFilename);
end

%% Determine frames to register
Config = load2PConfig(ImageFiles);

if strcmp(FrameIndex, 'new')
    if exist('MCdata', 'var')
        frameCounter = 0;
        FrameIndex = [];
        for findex = 1:numel(Config)
            FrameIndex = [FrameIndex, frameCounter+find(isnan(MCdata(findex).dpx(1, :, 1)))];
            frameCounter = frameCounter + Config(findex).Frames;
        end
    else
        FrameIndex = [1 inf];
    end
end

if FrameIndex(end) == inf
    FrameIndex = [FrameIndex(1:end-1), FrameIndex(end-1)+1:sum([Config(:).Frames])];
end

totalFrames = numel(FrameIndex);

%% Determine number batches
if LoadImgs
    numFiles = numel(ImageFiles);
    [Images, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', 2, 'Double');
    sizeFrame = whos('Images');
    sizeFrame = sizeFrame.bytes;
    if ispc
        mem = memory;
        nFramesPerLoad = max(1, floor(portionOfMemory*mem.MaxPossibleArrayBytes/sizeFrame));
    else
        nFramesPerLoad = max(1, floor(portionOfMemory*sizeRAM/sizeFrame));
    end
    numBatches = ceil(totalFrames/nFramesPerLoad);
else
    numFiles = 1;
    ImageFiles = {''};
    loadObj.Height = size(Images, 1);
    loadObj.Width = size(Images, 2);
    loadObj.Depth = size(Images, 3);
    loadObj.Channels = size(Images, 4);
    totalFrames = size(Images, 5);
    nFramesPerLoad = totalFrames;
    numBatches = 1;
end
fprintf('Performing fullDoLucasKanade on %d file(s) with %d frames...\n', numFiles, totalFrames);
fprintf('\t%s\n', ImageFiles{:}); % display each filename
fprintf('Calculation requires %d batches with %d frames per batch...\n', numBatches, nFramesPerLoad);

%% Initialize MCdata struct
if ~exist('MCdata', 'var')
    for findex = 1:numFiles
        MCdata(findex).type = 'doLucasKanade';
        MCdata(findex).date = datestr(now);
        MCdata(findex).FullFilename = ImageFiles{findex};
        MCdata(findex).Channel2AlignFrom = Channel2AlignFrom;
        MCdata(findex).Parameters = Parameters;
        MCdata(findex).dpx = nan(Parameters.Nbasis + 1, Config(findex).Frames, Config(findex).Depth);
        MCdata(findex).dpy = nan(Parameters.Nbasis + 1, Config(findex).Frames, Config(findex).Depth);
    end
else
    for findex = 1:numFiles
        MCdata(findex).type = 'doLucasKanade';
        MCdata(findex).date = datestr(now);
        MCdata(findex).FullFilename = ImageFiles{findex};
        MCdata(findex).Channel2AlignFrom = Channel2AlignFrom;
        MCdata(findex).Parameters = Parameters;
    end
end

%% Compute template
if ~Template
    if LoadImgs
        Images = load2P(ImageFiles, 'Type', 'Direct', 'Frames', 2:min(totalFrames, numFramesInitialTemplate+1), 'Double');
    end
    Template = mean(Images(:,:,:,Channel2AlignFrom,1:min(totalFrames, numFramesInitialTemplate)),5);
end

%% Initialize second pass parameters
if numFramesSecondPass
    MovingAvgFilt = 1/numFramesSecondPass*ones(1,numFramesSecondPass);
    AlignedFrames = zeros(loadObj.Height, loadObj.Width, loadObj.Depth, nFramesPerLoad);
    LastIter = [];
    Parameters2 = Parameters;
    Parameters2.deltacorr = Parameters2.deltacorr/2;
    Parameters2.niter = floor(Parameters2.niter/2);
end

%% Register images

if ischar(MCImgsFilename)
    fid = fopen(MCImgsFilename);
end

% Load each frames in batches and correct for motion
for bindex = 1:numBatches
    fprintf('\tbatch %d:', bindex);
    
    % Determine frames
    FirstFrame = (bindex-1)*nFramesPerLoad + 1;
    LastFrame = min(FirstFrame+nFramesPerLoad-1, totalFrames);
    if LoadImgs
        [Images, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', FrameIndex(FirstFrame):FrameIndex(LastFrame), 'Double');
    end
    numFrames = size(loadObj.FrameIndex, 1);
    
    % Register Each Frame to the Initial Template
    for currentFrame = 1:numFrames
        for currentDepth = 1:Config(findex).Depth
            [AlignedFrames(:,:,currentDepth,currentFrame),...
                MCdata(loadObj.FrameIndex(currentFrame,1)).dpx(:,loadObj.FrameIndex(currentFrame,2),currentDepth),...
                MCdata(loadObj.FrameIndex(currentFrame,1)).dpy(:,loadObj.FrameIndex(currentFrame,2),currentDepth)] = doLucasKanade_gpu(...
                    Template(:,:,currentDepth),...
                    Images(:,:,currentDepth,Channel2AlignFrom,currentFrame),...
                    zeros(Parameters.Nbasis+1, 1),...
                    zeros(Parameters.Nbasis+1, 1),...
                    Parameters);
%                     repmat(MCdata(loadObj.FrameIndex(currentFrame,1)).dpx(Parameters.Nbasis + 1,loadObj.FrameIndex(max(currentFrame-1,1),2),currentDepth), Parameters.Nbasis + 1, 1),... % seed with last frame's final segment's location as starting point
%                     repmat(MCdata(loadObj.FrameIndex(currentFrame,1)).dpy(Parameters.Nbasis + 1,loadObj.FrameIndex(max(currentFrame-1,1),2),currentDepth), Parameters.Nbasis + 1, 1),...
        end
    end
    fprintf('\tRegistered to template.');
    
    % Register Each Frame to a Local Average
    if numFramesSecondPass
        [AlignedFrames, LastIter] = filter(MovingAvgFilt, 1, AlignedFrames, LastIter, 4); % compute Moving-Average of Data
        if FirstFrame < numFramesSecondPass+1
            firstFrame = numFramesSecondPass+1;
        else
            firstFrame = 1;
        end
        for currentFrame = firstFrame:numFrames
            for currentDepth = 1:Config(findex).Depth
                [AlignedFrames(:,:,currentDepth,currentFrame),...
                    MCdata(loadObj.FrameIndex(currentFrame,1)).dpx(:,loadObj.FrameIndex(currentFrame,2),currentDepth),...
                    MCdata(loadObj.FrameIndex(currentFrame,1)).dpy(:,loadObj.FrameIndex(currentFrame,2),currentDepth)] = doLucasKanade_gpu(...
                        AlignedFrames(:,:,currentDepth,currentFrame),...
                        Images(:,:,currentDepth,Channel2AlignFrom,currentFrame),...
                        MCdata(loadObj.FrameIndex(currentFrame,1)).dpx(:,loadObj.FrameIndex(currentFrame,2),currentDepth),... % seed with results from first pass
                        MCdata(loadObj.FrameIndex(currentFrame,1)).dpy(:,loadObj.FrameIndex(currentFrame,2),currentDepth),...
                        Parameters2);
            end
        end
        fprintf('\tRegistered to local average.');
    end
    
    % Save motion corrected images
    if ischar(MCImgsFilename)
        for currentFrame = 1:numFrames
            for currentDepth = 1:Config(findex).Depth
                for c = 1:Config(findex).Channels
                    if c == Channel2AlignFrom % save already computed image
                        img = reshape(AlignedFrames(:,:,currentDepth,currentFrame), Config(findex).Height*Config(findex).Width, 1); % reshape back to a vector for saving
                    else % apply computed transformation to current channel
                        if findex==1 && bindex==1
                            [currentframe, B, xi, yi] = applyDoLucasKanade_gpu(...
                                Images(:,:,currentDepth,c,currentFrame),...
                                MCdata(loadObj.FrameIndex(currentFrame,1)).dpx(:,loadObj.FrameIndex(currentFrame,2),currentDepth),...
                                MCdata(loadObj.FrameIndex(currentFrame,1)).dpy(:,loadObj.FrameIndex(currentFrame,2),currentDepth));
                        else
                            currentframe = applyDoLucasKanade_gpu(...
                                Images(:,:,currentDepth,c,currentFrame),...
                                MCdata(loadObj.FrameIndex(currentFrame,1)).dpx(:,loadObj.FrameIndex(currentFrame,2),currentDepth),...
                                MCdata(loadObj.FrameIndex(currentFrame,1)).dpy(:,loadObj.FrameIndex(currentFrame,2),currentDepth), B, xi, yi);
                        end
                        img = reshape(currentframe, Config(findex).Height*Config(findex).Width, 1); % reshape back to a vector for saving
                    end
                    fwrite(fid, img, Config(1).Precision);
                end
            end
        end
        fprintf('\tFrames saved to file.');
    end
    
    % Save map
    if ischar(MCdataFilename)
        if bindex == 1 && ~exist(MCdataFilename, 'file')
            save(MCdataFilename, 'MCdata', '-mat');
        else
            save(MCdataFilename, 'MCdata', '-mat', '-append');
        end
        fprintf('\tTransform saved to file.');
    end
    
    fprintf('\n'); 
end %batch
fprintf('Registration complete\n');

% Notify of file saving
if ischar(MCImgsFilename)
    fprintf('Registered images saved to: %s\n', MCImgsFilename);
end
if ischar(MCdataFilename)
    fprintf('Registration map saved to: %s\n',MCdataFilename);
end


