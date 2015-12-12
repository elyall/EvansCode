function out = frameMean(Images, varargin)

Depth = 1;
Channel = 1;
FrameIndex = [1 inf];

% Memory settings
portionOfMemory = 0.08; % find 10% or less works best
sizeRAM = 32000000000; % amount of memory on your computer (UNIX-only)

directory = cd;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Depth'
                Depth = varargin{index+1};
                index = index + 2;
            case 'Channel'
                Channel = varargin{index+1};
                index = index + 2;
            case 'FrameIndex'
                FrameIndex = varargin{index+1};
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
    [f,p] = uigetfile('*.sbx;*.tif', 'Choose datasets', directory, 'MultiSelect', 'on');
    if isnumeric(f)
        return
    end
    Images = fullfile(p,f);
elseif ischar(Images)
    Images = {Images};
end


%% Load in Data and determine dimensions
if iscellstr(Images) % filename input
    loadType = true;
    ImageFile = Images;
    [Images, loadObj, Config] = load2P(ImageFile, 'Type', 'Direct', 'Frames', 2, 'Double', 'Depth', Depth, 'Channel', Channel);
    sizeFrame = whos('Images');
    sizeFrame = sizeFrame.bytes;
    if ispc
        mem = memory;
        numFramesPerLoad = max(1, floor(portionOfMemory*mem.MaxPossibleArrayBytes/sizeFrame));
    else
        numFramesPerLoad = max(1, floor(portionOfMemory*sizeRAM/sizeFrame));
    end
    numFrames = sum([Config(:).Frames]);
    dim = loadObj.size;
else % numeric array input
    loadType = false;
    dim = size(Images);
    numFrames = size(Images, 5);
    numFramesPerLoad = numFrames; % all frames already loaded
end


%% Determine frames
if FrameIndex(end) == inf
    FrameIndex = [FrameIndex(1:end-1),FrameIndex(1:end-1)+1:numFrames];
end
totalFrames = numel(FrameIndex);


%% Take mean

out = nan(totalFrames,1); % intialize output

% Load frames in batches and compute mean
for bindex = 1:numFramesPerLoad:totalFrames % direct loading only -> load frames in batches
    lastframe = min(bindex+numFramesPerLoad-1, totalFrames);
    currentFrames = FrameIndex(bindex:lastframe);
    
    % Load frames
    if loadType
        Images = load2P(ImageFile, 'Type', 'Direct', 'Frames', currentFrames, 'Double', 'Depth', Depth, 'Channel', Channel); %direct
    end
    
    % Compute mean
    out(bindex:lastframe) = mean(reshape(Images, dim(1)*dim(2), numel(currentFrames)));
    
end

