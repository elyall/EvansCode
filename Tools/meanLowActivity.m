function Mean = meanLowActivity(Images, AnalysisInfo, varargin)

NumExtraFrames = 8;
MotionCorrect = false;
% FrameIndex = [];
% TrialIndex = [1 inf];

Depth = 1;
Channel = 1;

% Memory settings
portionOfMemory = 0.08; % find 10% or less works best
sizeRAM = 32000000000; % amount of memory on your computer (UNIX-only)

directory = cd;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
%             case {'Trials','trials','TrialIndex'}
%                 TrialIndex = varargin{index+1};
%                 index = index + 2;
            case 'NumExtraFrames'
                NumExtraFrames = varargin{index+1};
                index = index + 2;
            case 'MotionCorrect'
                MotionCorrect = varargin{index+1};
                index = index + 2;
%             case 'FrameIndex'
%                 FrameIndex = varargin{index+1};
%                 index = index + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end


%% Load in images
if iscellstr(Images) || ischar(Images) % filename input
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


%% Motion correction
if ischar(MotionCorrect)
    load(MotionCorrect, 'MCdata', '-mat');
    MotionCorrect = true;
elseif isstruct(MotionCorrect)
    MCdata = MotionCorrect;
    MotionCorrect = true;
elseif MotionCorrect == true && ~exist('MCdata', 'var')
    load(ExperimentFile, 'MCdata', '-mat');
end


%% Determine frames with low activity
% IndivMean = frameMean(Images);
% Index = IndivMean <= prctile(IndivMean, 30);
% % Index(~ismember(Index,FrameIndex)) = []; % remove unwanted frames


%% Determine frames that should have low activity
if ischar(AnalysisInfo)
    load(AnalysisInfo, 'AnalysisInfo', '-mat');
end

Index = setdiff(2:numFrames,cell2mat(arrayfun(@colon, AnalysisInfo.ExpStimFrames(:,1), AnalysisInfo.ExpStimFrames(:,2)+NumExtraFrames, 'Uniform', false)'));
% Index(~ismember(Index,FrameIndex)) = []; % remove unwanted frames


%% Load in frames and take average
totalFrames = numel(Index);

Mean = zeros(dim(1:2));

% Load frames in batches and compute mean
for bindex = 1:numFramesPerLoad:totalFrames % direct loading only -> load frames in batches
    lastframe = min(bindex+numFramesPerLoad-1, totalFrames);
    currentFrames = Index(bindex:lastframe);
    
    % Load frames
    if loadType
        [Images, loadObj] = load2P(ImageFile, 'Type', 'Direct', 'Frames', currentFrames, 'Double', 'Depth', Depth, 'Channel', Channel); %direct
        if MotionCorrect
            Images = applyMotionCorrection(Images, MCdata, loadObj);
        end
    end
    
    % Compute mean
    Mean = Mean + sum(bsxfun(@rdivide, Images, totalFrames),5);
    
end


