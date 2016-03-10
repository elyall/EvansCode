function timestamps = loadAVITimeStamps(VideoFile, varargin)
%(f,1)= second (0 to 127)
%(f,2)= cycle_count (0 to 7999)
%(f,3)= cycle_offset (0 to x -> x depends on configuration)

offset = 0;
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
            case 'offset'
                offset = varargin{index+1};
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

JUNK         = [74 85 78 75];			% <JUNK>
MOVI         = [109 111 118 105];		% <movi>
ValidFrameID = [48 48 100 98];			% <00db>
STRF         = [115 116 114 102];		% <strf>

%% Load in file header information

% Load in file header
fid = fopen(VideoFile, 'r');
if fid < 3; error('File ''%s'' not found',VideoFile); end
Header = fread(fid, 100000, 'uint8')';
fclose(fid);

% Perform sanity checks
if ~all(Header(1:4)==uint8('RIFF'))										%check for: 'RIFF'
    error('%s not a valid RIFF file.',VideoFile);
elseif ~all(Header(9:11)==uint8('AVI'));									%check for: 'AVI '
    error('%s not a valid AVI file.',VideoFile);
end

% Determine frame information
o = findstr(Header,uint8('avih')); o=o(1);
% time_per_frame = sum(Header(o+8:o+11).*[1,256,256^2,256^3]);
numFrames      = sum(Header(o+24:o+27).*[1,256,256^2,256^3]);
Width          = sum(Header(o+40:o+43).*[1,256,256^2,256^3]);
Height         = sum(Header(o+44:o+47).*[1,256,256^2,256^3]);

% Determine bit depth
% o = findstr(Header,uint8('strf')); o=o(1);
% color_depth=Header(o+22)+Header(o+23)*256;
% switch color_depth
%    case 8
%       bytes=1;
%    case 16
%       bytes=2;
%    case 24 
%       bytes=3;
%    otherwise bytes=4;
% end

% Determine frame identifier info (based on first frame found)
o = findstr(Header,uint8('movi')); o=o(1);
frame_ID=Header(o+4:o+7);


%% Determine frames to read in
if Frames(end) == inf
    Frames = [Frames(1:end-1),Frames(end-1)+1:numFrames];
end
diffFrames = diff([0,Frames])'-1;                               % number of frames between each subsequent frame to be loaded
seekOperations = find(diffFrames);                              % find jumps
seekOperations = [seekOperations, diffFrames(seekOperations)];  % 
numFramesToLoad = numel(Frames);


%% Load in requested frames
fid = fopen(VideoFile, 'r');
if fid < 3; error('File ''%s'' not found',VideoFile); end
fseek(fid,o+3,-1); % move to beginning of first frame

timestamps = nan(numFramesToLoad,3); %initialize output

% Load in information from requested frames
for findex = 1:numFramesToLoad
    
    if isempty(seekOperations) || ~ismember(seekOperations(:,1), findex)
        numSkip = 0;
    else
        numSkip = seekOperations(seekOperations(:,1)==findex,2);
    end
   
    while 1
        
        frame_header = fread(fid, 8, 'uint8')';                 % load in current frame's header
        f_length = sum(frame_header(5:8).*[1,256,256^2,256^3]); % determine frame's length
        
        if isequal(frame_header(1:4),frame_ID) && numSkip       % found valid frame, but data not requested
            disp('skipping unwanted frame...');
            fseek(fid, f_length, 0);                            % skip frame
            numSkip = numSkip - 1;
            
        elseif f_length==0;                                     % found empty frame -> move on
            disp('skipping empty frame...');
            
        elseif isequal(frame_header(1:4),frame_ID)              % found valid frame
            break
            
        elseif ~isequal(frame_header(1:4),frame_ID)             % found non-movie frame
            %disp('skipping non-movie frame...');
            %fseek(fid, f_length, 0);                            % skip frame
            %f_length
            
        elseif isequal(frame_header(1:4),uint8('JUNK'))         % found JUNK frame
            disp('skipping JUNK frame...');
            fseek(fid, f_length, 0);                            % skip junk
            
        end

    end
    
    % move to last column of frame
    fseek(fid, Width*(Height-1), 0);
    
    % timestamp information is not at front of frame
    if offset
        fseek(fid, offset, 0);
    end
    
    % load in timestamp information
    timestamps(findex, 1) = fread(fid, 1, 'ubit7');     % second (0 to 127)
    timestamps(findex, 2) = fread(fid, 1, 'ubit13');    % cycle_count (0 to 7999)
    timestamps(findex, 3) = fread(fid, 1, 'ubit12');    % cycle_offset (0 to x -> x depends on configuration)
    
    fseek(fid, Width-4, 0);                          % move to end of frame
    
end
fclose(fid); % close file

    
