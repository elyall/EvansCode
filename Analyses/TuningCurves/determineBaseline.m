function Baseline = determineBaseline(Images,Prctile,varargin)

Filter = [];
% Filter = fspecial('gaussian', 5, 1);
Frames = 2:2000;    % vector of frame indices to load (loading only)
Depth = [1,inf];    % vector of depth indices to load (loading only)
Channel = 1;        % vector of channel indices to load (loading only)
MCdata = [];        % MCdata object for motion-correcting frames (loading only)

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Filter'
                Filter = varargin{index+1};
                index = index + 2;
            case 'Frames'
                Frames = varargin{index+1};
                index = index + 2;
            case 'Depth'
                Depth = varargin{index+1};
                index = index + 2;
            case 'Channel'
                Channel = varargin{index+1};
                index = index + 2;
            case 'MCdata'
                MCdata = varargin{index+1};
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
    [Images,directory] = uigetfile({'*.sbx;*.tiff'}, 'Select image files:', directory,'MultiSelect','on');
    if isnumeric(Images)
        return
    end
    Images = fullfile(directory,Images);
end

if ~exist('Prctile','var') || isempty(Prctile)
    Prctile = 30;
end

if ischar(MCdata)
    load(MCdata,'MCdata','-mat');
end

fprintf('Determining baseline...');

%% Load images
if iscellstr(Images) || ischar(Images)
    Config = load2PConfig(Images);
    if Frames(end) > sum(Config(:).Frames)
        Frames = Frames(1):sum(Config(:).Frames);
    end
    [Images,loadObj] = load2P(Images,'Frames',Frames,'Channel',Channel,'Depth',Depth,'Double');
    if ~isempty(MCdata)
        Images = applyMotionCorrection(Images,MCdata,loadObj);
    end
end

%% Blur images
if ~isempty(Filter)
    Images = imfilter(Images,Filter);
end

%% Compute baseline image
Baseline = prctile(Images,Prctile,5);

fprintf('\tComplete\n');
