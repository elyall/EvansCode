function config = parseSbxHeader(File)
% File can be a .mat info file, a .sbx file, a directory to prompt from, or
% an empty matrix to initiate file selection from the default directory

% config.version = 1;

%% Check input arguments
narginchk(0,1);
if ~exist('File', 'var') || isempty(File)
    directory = loadCanalSettings('DataDirectory');
    [File, p] = uigetfile({'.sbx;.mat'}, 'Choose sbx file', directory);
    if isnumeric(File)
        return
    end
    File = fullfile(p, File);
elseif isdir(File)
    [File, p] = uigetfile({'.sbx;.mat'}, 'Choose sbx file', File);
    if isnumeric(File)
        return
    end
    File = fullfile(p, File);
end


%% Load in header
[~, ~, e] = fileparts(File);
switch e
    case '.mat'
        ConfigFile = File;
        temp = sbxIdentifyFiles(File);
        SbxFile = temp{1};
    case '.sbx'
        SbxFile = File; % assumes sbx file to have same name and be located on same path
        temp = sbxIdentifyFiles(File);
        ConfigFile = temp{1};
end


%% Set identifying info
config.type = 'sbx';
config.FullFilename = SbxFile;
[~, config.Filename, ~] = fileparts(config.FullFilename);
config.Tag = getTag(config.Filename);


%% Identify header information from file

% Load header
load(ConfigFile, 'info');

% Update header
if ~isfield(info, 'Width')
    info = updateInfoFile(SbxFile, ConfigFile);
end

% Save header
config.header = {info};

% Save important information
config.Height = info.Height;
config.Width = info.Width;

% Determine # of channels and # of frames
config.Channels = info.numChannels;
config.Frames = info.numFrames;

% Determine magnification
config.ZoomFactor = info.config.magnification;

%% DEFAULTS
% config.Processing = {};
% config.info = [];
% config.MotionCorrected = false;
config.FrameRate = 15.5; % current default
config.Depth = 1; % current default
config.ZStepSize = 0; % current default
config.Precision = 'uint16'; % default
config.DimensionOrder = {'Channels','Width','Height','Frames','Depth'}; % default
config.Colors = {'green', 'red'};
config.size = sizeDimensions(config);
