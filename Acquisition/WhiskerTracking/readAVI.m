function [data, timestamps, framecounter] = readAVI(VideoFile, varargin)

Frames = [1 inf];
timestamps = 1;
framecounter = 4; 

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


%% Load in requested frames

% Initialize output
if vidObj.BitsPerPixel == 8
    Precision = 'uint8';
end
switch vidObj.VideoFormat
    case 'Grayscale'
        numChannels = 1;
end
data = zeros(vidObj.Height,vidObj.Width,numChannels,numFramesToLoad,Precision);

% Load in frames
for findex = 1:numFramesToLoad
    if vidObj.CurrentTime ~= Frames(findex)     % video object not currently in right place
        vidObj.CurrentTime = Frames(findex);    % move to correct time
    end
    data(:,:,:,findex) = readFrame(vidObj);     % read in current frame
end


%% Determine timestamps
if timestamps
    temp = double(squeeze([data(1,timestamps:timestamps+3,1,:)])');
    
    % Convert 4xuint8 to [uint7, uint13, uint12] (uint12 is uint8 for USB)
    temp = de2bi(temp',8,'left-msb');
    temp = reshape(temp',8*4,numFramesToLoad)';
    temp = [bi2de(temp(:,1:7),2,'left-msb'), bi2de(temp(:,8:20),2,'left-msb'), bi2de(temp(:,21:28),2,'left-msb')];
    
    timestamps = [zeros(1,3);diff(temp,[],1)]; %convert to difference from cyclical
    
    % Fix jumps in cycle
    bits = [7,13,8];
    for cindex = 1:3
        jumps = find(timestamps(:,cindex)<0);
        for lindex = 1:numel(jumps)
            timestamps(jumps(lindex),cindex) = temp(jumps(lindex),cindex)+2^bits(cindex)-temp(jumps(lindex)-1,cindex);
        end
    end
end


%% Determine frame counter
if framecounter
    temp = double(squeeze([data(1,framecounter:framecounter+3,1,:)])');
    
    framecounter = [zeros(1,4);diff(temp,[],1)]; %convert to difference from cyclical
    
    % Fix jumps in cycle
    for cindex = 1:4
        jumps = find(framecounter(:,cindex)<0);
        for lindex = 1:numel(jumps)
            framecounter(jumps(lindex),cindex) = temp(jumps(lindex),cindex)+256-temp(jumps(lindex)-1,cindex);
        end
    end
end

