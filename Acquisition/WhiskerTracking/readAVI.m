function data = readAVI(VideoFile, varargin)

Frames = [1 inf];

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Frames','frames','Frame','frame'}
                Frames = varargin{index+1};
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

if ~exist('VideoFile','var') || isempty(VideoFile)
    [VideoFile,p] = uigetfile({'*.avi'},'Select video file to load',directory);
    if isnumeric(VideoFile)
        return
    end
    VideoFile = fullfile(p,VideoFile);
end

vidObj = VideoReader(VideoFile);


%% Determine frames to read in
numFrames = vidObj.Duration*vidObj.FrameRate;
if Frames(end) == inf
    Frames = [Frames(1:end-1),Frames(end-1)+1:numFrames];
end
Frames = (Frames-1)/vidObj.FrameRate; % convert to time indexing
numFramesToLoad = numel(Frames);

% % Determine number and location of seek operations due to non-contiguous frame requests (jumps)
% seekoperations = find(diff(Frames)~=1); %find any jumps within the frames to read out (jumps requiring seeking)
% if isempty(seekoperations) %no jumps
%     numframesperread = numFramesToLoad; %all frames will be read in one read
%     seekoperations = 1; %only one seek operation to first frame of FrameIndex
% else
%     numframesperread = diff([0,seekoperations,numFramesToLoad]); %multiple reads required with various numbers of frames per read
%     seekoperations = [1,seekoperations+1]; %indexes the first frame of each read within FrameIndex
% end


%% Initialize output
if vidObj.BitsPerPixel == 8
    Precision = 'uint8';
end
switch vidObj.VideoFormat
    case 'Grayscale'
        numChannels = 1;
end

data = zeros(vidObj.Height,vidObj.Width,numChannels,numFramesToLoad,Precision);


%% Load in requested frames
for findex = 1:numFramesToLoad
    if vidObj.CurrentTime ~= Frames(findex)     % video object not currently in right place
        vidObj.CurrentTime = Frames(findex);    % move to correct time
    end
    data(:,:,:,findex) = readFrame(vidObj);     % read in current frame
end


