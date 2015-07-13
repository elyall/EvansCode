function [MCdata, Template] = fullDoLucasKanade1(Images, Template, varargin)

% Default settings
Parameters.Nbasis = 16;
Parameters.niter = 25;
Parameters.damping = 1;
Parameters.deltacorr = .0005;

Channel2AlignFrom = 1;
numFramesInitialTemplate = 500;
numFramesSecondPass = 20;
MCdataFilename = false;
MCImgsFilename = false;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    switch varargin{index}
        case 'SaveAlignmentTo'
            MCdataFilename = varargin{index+1};
            index = index + 2;
        case 'SaveImagesTo'
            MCImgsFilename = varargin{index+1};
            index = index + 1;
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
    [f,p] = uigetfile({'*.mat'}, 'Save alignment data to:', directory);
    if isnumeric(f)
        return
    end
    MCdataFilename = fullfile(p, f);
end
if MCImgsFilename == true
    directory = CanalSettings('DataDirectory');
    [f,p] = uiputfile({'*.sbx'}, 'Save registered images as:', directory);
    if isnumeric(f)
        return
    end
    MCImgsFilename = fullfile(p, f);
end

%% Determine # of frames and initialize struct
if LoadImgs
    [Images, config] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', 2, 'Double');
    mem = memory;
    sizeFrame = whos('Images');
    sizeFrame = sizeFrame.bytes;
    nFramesPerLoad = max(1, floor(0.08*mem.MaxPossibleArrayBytes/sizeFrame));
    fprintf('Aligning %d frames. \tFilename: %s', config.Frames, ImageFiles{1});
else
    config.Height = size(Images, 1);
    config.Width = size(Images, 2);
    config.Depth = size(Images, 3);
    config.Channels = size(Images, 4);
    nFramesPerLoad = size(Images, 5);
    config.Frames = nFramesPerLoad;
    ImageFiles = ''; % for waitbar
    fprintf('Aligning %d frames', config.Frames);
end

if ~exist('MCdata', 'var')
    MCdata.type = 'doLucasKanade';
    MCdata.dpx = zeros(Parameters.Nbasis + 1, config.Frames, config.Depth);
    MCdata.dpy = zeros(Parameters.Nbasis + 1, config.Frames, config.Depth);
end

%% Compute template and filter info
if ~Template
    if LoadImgs
        Images = load2P(ImageFiles, 'Type', 'Direct', 'Frames', 1:min(config.Frames, numFramesInitialTemplate), 'Double');
    end
    Template = mean(Images(:,:,:,Channel2AlignFrom,1:min(config.Frames, numFramesInitialTemplate)),5);
end

% Second pass info
if numFramesSecondPass
    MovingAvgFilt = 1/numFramesSecondPass*ones(1,numFramesSecondPass);
    AlignedFrames = zeros(config.Height, config.Width, config.Depth, nFramesPerLoad);
    LastIter = [];
    Parameters2 = Parameters;
    Parameters2.deltacorr = Parameters2.deltacorr/2;
    Parameters2.niter = floor(Parameters2.niter/2);
end

%% Register images
nBatches = ceil(config.Frames/nFramesPerLoad);
BatchIndex = 1;
wb=waitbar(1/config.Frames, sprintf('dLK: batch %d of %d - first pass',BatchIndex,nBatches));
for f = 1:nFramesPerLoad:config.Frames
    
    % Determine frames in batch
    lastframe = min(f+nFramesPerLoad-1, config.Frames);
    if LoadImgs
        Images = load2P(ImageFiles, 'Type', 'Direct', 'Frames', f:lastframe, 'Double');
    end
    FrameIndices = f:lastframe;
    
    % Register Each Frame to the Initial Template
    waitbar(f/config.Frames, wb, sprintf('dLK: batch %d of %d - first pass',BatchIndex,nBatches));
    for currentFrame = 1:lastframe-f+1
        for currentDepth = 1:config.Depth
               [AlignedFrames(:,:,currentDepth,currentFrame),...
                MCdata.dpx(:,FrameIndices(currentFrame),currentDepth),...
                MCdata.dpy(:,FrameIndices(currentFrame),currentDepth)] = doLucasKanade_gpu(Template,...
                Images(:,:,currentDepth,Channel2AlignFrom,currentFrame),...
                repmat(MCdata.dpx(Parameters.Nbasis + 1,FrameIndices(max(currentFrame-1,1)),currentDepth), Parameters.Nbasis + 1, 1),... % use last frame's final segment's location as starting point
                repmat(MCdata.dpy(Parameters.Nbasis + 1,FrameIndices(max(currentFrame-1,1)),currentDepth), Parameters.Nbasis + 1, 1),...
                Parameters);    
        end
    end
    
    % Register Each Frame to a Local Average
    if numFramesSecondPass
        waitbar(f/config.Frames, wb, sprintf('dLK: batch %d of %d - second pass',BatchIndex,nBatches));
        [AlignedFrames, LastIter] = filter(MovingAvgFilt, 1, AlignedFrames, LastIter, 4); % compute Moving-Average of Data
        if f < numFramesSecondPass+1
            firstFrame = numFramesSecondPass+1;
        else
            firstFrame = 1;
        end
        for currentFrame = firstFrame:lastframe-f+1
            for currentDepth = 1:config.Depth
                [AlignedFrames(:,:,currentDepth,currentFrame),...
                    MCdata.dpx(:,FrameIndices(currentFrame),currentDepth),...
                    MCdata.dpy(:,FrameIndices(currentFrame),currentDepth)] = doLucasKanade_gpu(AlignedFrames(:,:,currentDepth,currentFrame),...
                    Images(:,:,currentDepth,Channel2AlignFrom,currentFrame),...
                    MCdata.dpx(:,FrameIndices(currentFrame),currentDepth),...
                    MCdata.dpy(:,FrameIndices(currentFrame),currentDepth),...
                    Parameters2);
            end
        end
    end
    
    % Save motion corrected images
    if MCImgsFilename
        waitbar(f/config.Frames, wb, sprintf('dLK: batch %d of %d - saving batch',BatchIndex,nBatches));
        for currentFrame = 1:lastframe-f+1
            for currentDepth = 1:config.Depth
                for c = 1:config.Channels
                    if c == Channel2AlignFrom
                        % Save already computed image
                        img = reshape(AlignedFrames(:,:,currentDepth,currentFrame), config.Height*config.Width, 1); % reshape back to a vector for saving
                        fwrite(fid,img,config.Precision);
                    else
                        % Align frame from current channel
                        if f == 1
                            [currentframe, B, xi, yi] = applyDoLucasKanade_gpu(Images(:,:,currentDepth,c,currentFrame), MCdata.dpx(:,FrameIndices(currentFrame),currentDepth), MCdata.dpy(:,FrameIndices(currentFrame),currentDepth));
                        else
                            currentframe = applyDoLucasKanade_gpu(Images(:,:,currentDepth,c,currentFrame), MCdata.dpx(:,FrameIndices(currentFrame),currentDepth), MCdata.dpy(:,FrameIndices(currentFrame),currentDepth), B, xi, yi);
                        end
                        img = reshape(currentframe,config.Height*config.Width,1); % reshape back to a vector for saving
                        fwrite(fid,img,config.Precision);
                    end
                end
            end
        end
    end
    
    BatchIndex = BatchIndex + 1;
end
close(wb);

% Notify of file saving
if MCImgsFilename
    fprintf('Registered images saved to: %s\n', ImgsFilename);
end

% Format output
MCdata.Parameters = Parameters;
MCdata.Channel2AlignFrom = Channel2AlignFrom;
% Save map
if MCdataFilename
    fprintf('Saving registration map to file: %s...',MCdataFilename);
    save(MCdataFilename, 'MCdata', '-mat', '-append');
    fprintf('\tComplete\n');
end


