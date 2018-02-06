function Images = avgImgsN(Images,N,MCdata)
% AVGIMGSN produces mean image by averaging largest or smallest N pixels

direction = 'descend';
FrameIndex = [1,inf];
Depth = 3;
saveOut = true;
saveFile = 'temp3.tif';

% Memory settings
portionOfMemory = 0.08;     % find 10% or less works best
sizeRAM = 32000000000;      % amount of memory on your computer (UNIX-only)


%% Parse input arguments & load data
if ~exist('Images','var')
    Images = [];
end
if ischar(Images) || iscellstr(Images) || isempty(Images)
    loadImgs = true;
    ImageFiles = Images;
    [Images,~,Config] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', 2, 'Double');
    numFrames = Config.Frames;
    sizeFrame = whos('Images');
    sizeFrame = sizeFrame.bytes;
    if ispc
        mem = memory;
        numFramesPerLoad = max(1, floor(portionOfMemory*mem.MaxPossibleArrayBytes/sizeFrame));
    else
        numFramesPerLoad = max(1, floor(portionOfMemory*sizeRAM/sizeFrame));
    end
%     numFramesPerLoad = numFramesPerLoad-rem(numFramesPerLoad,Config.Depth);
    numFramesPerLoad = numFramesPerLoad*Config.Depth;
else
    loadImgs = false;
    numFrames = size(Images,5);
end

if ~exist('N','var') || isempty(N)
    N = 1000;
end

if ~exist('MCdata','var')
    MCdata = [];
elseif ischar(MCdata)
    load(MCdata,'MCdata','-mat');
end

if FrameIndex(end)==inf
    FrameIndex = cat(2, FrameIndex(1:end-1), FrameIndex(end-1)+1:numFrames);
end


%% Compute average
if loadImgs
    Images = [];
    for bindex = 1:numFramesPerLoad:numFrames                   % load frames in batches
        fprintf('On %d of %d frames\n',bindex,numFrames);
        lastframe = min(bindex+numFramesPerLoad-1, numFrames);  % determine if reached end
        currentFrames = FrameIndex(bindex:lastframe);           % determine frames to load
        [frames, loadObj] = load2P(ImageFiles, 'Type', 'Direct', 'Frames', currentFrames, 'Depth', Depth); % load current frames
        if ~isempty(MCdata)
            frames = applyMotionCorrection(frames, MCdata, loadObj); % apply motion correction
        end
        Images = sort(cat(5,Images,frames),5,direction); % sort frames
        Images = Images(:,:,:,:,1:min(N,end));           % keep only # of frames desired
    end
    
else
    if ~isempty(MCdata)
        Images = applyMotionCorrection(Images, MCdata); % apply motion correction
    end
    Images = Images(:,:,:,:,FrameIndex); % keep only requested frames
    Images = sort(Images,5,direction);   % sort frames
    Images = Images(:,:,:,:,1:N);        % keep only # of frames desired
end

Images = mean(Images,5); % compute average


%% Save output
if saveOut
    save2P(saveFile,Images);
end