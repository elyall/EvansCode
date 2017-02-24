function Avg = avgStack(Images,numFramesPerDepth,varargin)

MCdata = [];
numDepths = [];
saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'MCdata'
                MCdata = varargin{index+1};
                index = index + 2;
            case 'numDepths'
                numDepths = varargin{index+1};
                index = index + 2;
            case {'Save','save'}
                saveOut = true;
                index = index + 1;
            case {'saveFile','SaveFile'}
                saveFile = varargin{index+1};
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

if ~exist('Images','var') || isempty(Images)
    [Images,p] = uigetfile({'*.tif;*.Sbx'}, 'Select files containing maps:', directory);
    if isnumeric(Images)
        return
    end
    Images = fullfile(p,Images);
end

if ~exist('numFramesPerDepth','var') || isempty(numFramesPerDepth)
    numFramesPerDepth = inf;
end


%% Load data
if ischar(Images)
    ImagesFile = Images;
    Images = load2P(ImagesFile, 'Frames', [2,inf]);
    if saveOut && isempty(saveFile)
        saveFile = [ImagesFile(1:end-4),'_avg.tif'];
    end
end


%% Determine number of frames per depth & number of depths
[H,W,D,C,F] = size(Images);
if isinf(numFramesPerDepth)
    numFramesPerDepth = F;
end
if isempty(numDepths);
    numDepths = floor(F/numFramesPerDepth(end));
end
if numel(numFramesPerDepth)<numDepths
    try
        numFramesPerDepth = [numFramesPerDepth;repmat(numFramesPerDepth(end),numDepths-numel(numFramesPerDepth),1)];
    catch
        numFramesPerDepth = [numFramesPerDepth';repmat(numFramesPerDepth(end),numDepths-numel(numFramesPerDepth),1)];
    end
end


%% Motion correct frames
if ~isempty(MCdata)
    Images = applyMotionCorrection(Images,MCdata);
end


%% Average across frames
Avg = nan(H,W,numDepths);
for dindex = 1:numDepths
    if dindex ~= numDepths
        Avg(:,:,dindex) = mean(Images(:,:,1,1,numFramesPerDepth*(dindex-1)+1:numFramesPerDepth*dindex),5);
    else
        Avg(:,:,dindex) = mean(Images(:,:,1,1,numFramesPerDepth*(dindex-1)+1:end),5);
    end
end


%% Save output
if saveOut
    save2P(saveFile, Avg, 'CLim', true);
end

