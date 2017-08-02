function [Filename, CLim, indMap, Dim, Map] = vidStim(Filename, Images, varargin)
%vid Saves images to a video file
%   FILENAME = vidStim() prompts user to select a FILENAME to save to
%   and select a .exp or .align file to load data from.
%
%   FILENAME = vidStim(FILENAME, IMAGES) saves IMAGES to the video file
%   FILENAME. IMAGES can be a cell array of filenames to load data from, or
%   a cell array of cell arrays each containing a 3-dimensional matrix of
%   image data [H x W x F].
%


% Default parameters that can be adjusted

% Stim specific
StimIndex = [2,inf];        % vector specifying indices of stimuli to save
ControlIndex = false;       % index of control stimulus within cell array after accounting stims removed due to StimIndex
VarToLoad = 'AvgTrialdFoF'; % variable name to load from IMAGES if IMAGES is a char or cell array of strings
StimFrameIndex = [];        % numStim x 2 array specifying frame indices of first and last frame within each stimulus

% Overlays & colorbar
Text = {};                  % cell array of text to overlay
TextIndex = [];             % array of length numFrames indexing which text to display on each frame
FontSize = 30;              % scalar specifying size of text
Color = [1,1,1];            % color of overlays
showColorBar = false;       % whether to display the colorbar
CbLabel = 'dF/F'; 
CbFontSize = 20; 

% Specify data
% MCdata = {[]};              % cell array of MCdata structures
Crop = false;               % false or [numFiles x 4] vector specifying number of pixels to remove from edges (see: crop)
borderLims = false;         % false or [numFiles x 4] vector specifying number of pixels along edges to set to zero, inclusive (top, bottom, left, right) (faster than crop)
Afilt = false;               % 2D filter to apply to each frame
% Afilt = fspecial('gaussian',5,1);               % 2D filter to apply to each frame
Tfilt = false;              % 1D filter to apply across time
% Tfilt = 1/100*ones(1,100);              % 1D filter to apply across time

% Merging properties
Maps = [];                  % cell array or array of imref2d objects specifying location of each input dataset
MergeType = 'mean';         % 'mean' or 'blend' specifying how to merge the datasets ('blend' takes a long time)
% Only used if want to skip repeating long blend process when making
% multiple videos (see: mapFoVs)
indMap = [];                % H x W x numFiles indexing array (see: mapFoVs)
Dim = [];                   % numFiles x 4 specifying the locations and size of each dataset (see: mapFoVs)
Map = [];                   % imref2d object of the composite data (see: mapFoVs)

% Display properties
CMap = 'parula';            % Nx3 colormap, or string specifying the colormap of the video
CLim = [];                  % 1x2 vector specifying the color limits ([] uses max and min)
frameRate = 15.45*2;        % scalar specifying the frame rate of the output video
pixelSize = [1,1];          % 1x2 vector specifying the aspect ratio of the Height and Width
outputSize = [];            % 1x2 vector specifying the desired Height and Width of the output video in pixels


% Placeholders
directory = cd; % default directory when prompting user to select a file

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case 'ControlIndex'
                ControlIndex = varargin{index+1};
                index = index + 2;
            case 'VarToLoad'
                VarToLoad = varargin{index+1};
                index = index + 2;
            case 'StimFrameIndex'
                StimFrameIndex = varargin{index+1};
                index = index + 2;
            case 'Text'
                Text = varargin{index+1};
                index = index + 2;
            case 'TextIndex'
                TextIndex = varargin{index+1};
                index = index + 2;
            case 'FontSize'
                FontSize = varargin{index+1};
                index = index + 2;
            case 'Color'
                Color = varargin{index+1};
                index = index + 2;
            case {'ColorBar','Colorbar','colorbar'}
                showColorBar = varargin{index+1};
                index = index + 2;
            case 'CbLabel'
                CbLabel = varargin{index+1};
                index = index + 2;
            case 'CbFontSize'
                CbFontSize = varargin{index+1};
                index = index + 2;
%             case 'MCdata'
%                 MCdata = varargin{index+1};
%                 index = index + 2;
            case 'Crop'
                Crop = varargin{index+1};
                index = index + 2;
            case 'borderLims'
                borderLims = varargin{index+1};
                index = index + 2;
            case {'Afilter','Afilt'}
                Afilt = varargin{index+1};
                index = index + 2;
            case {'Tfilter','Tfilt'}
                Tfilt = varargin{index+1};
                index = index + 2;
            case 'Maps'
                Maps = varargin{index+1};
                index = index + 2;
            case 'MergeType'
                MergeType = varargin{index+1};
                index = index + 2;
            case 'Map'
                indMap = varargin{index+1};
                Dim = varargin{index+2};
                Map = varargin{index+3};
                index = index + 4;
            case 'CMap'
                CMap = varargin{index+1};
                index = index + 2;
            case 'CLim'
                CLim = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
                index = index + 2;
            case 'pixelSize'
                pixelSize = varargin{index+1};
                index = index + 2;
            case 'outputSize'
                outputSize = varargin{index+1};
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

if ~exist('Filename', 'var') || isempty(Filename)
    [Filename, p] = uiputfile({'*.avi'},'Save file as',directory);
    if isnumeric(Filename)
        return
    end
    Filename = fullfile(p, Filename);
end

if ~exist('Images', 'var') || isempty(Images)
    [Images,p] = uigetfile({'*.exp;*.align'}, 'Select image files for first dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(Images)
        return
    end
    Images = fullfile(p, Images);
end
if ~iscell(Images)
    Images = {{Images}};
end
numFiles = numel(Images);
for findex = 1:numFiles
    if ~iscell(Images{findex})
        Images{findex} = {Images{findex}};
    end
end


%% Load and format images
if iscellstr(Images)
    ImageFiles = Images;
    Images = cell(numel(ImageFiles),1);
    for findex = 1:numel(ImageFiles)
        temp = load(ImageFiles{findex}, VarToLoad, '-mat');
        Images{findex} = temp.(VarToLoad);
    end
end
N = ndims(Images{1}{1}); % determine frame dimension
numFrames = cellfun(@(x) size(x,N), Images{1}); % number of frames in each cell

%% Keep and concatenate requested stimuli

% Determe stimuli to keep
if StimIndex(end) == inf
    StimIndex = [StimIndex(1:end-1), StimIndex(end-1)+1:numel(Images{1})];
end
numStims = numel(StimIndex);

% Concatenate requested stimuli
for findex = 1:numFiles
    Images{findex} = cat(N,Images{findex}{StimIndex});
end
numFrames(setdiff(1:numel(numFrames),StimIndex)) = []; % throw out removed stimuli


%% Format stim frame index

% Replicate frame index for all stimuli
if size(StimFrameIndex,1) == 1 && numel(StimIndex) > 1
    StimFrameIndex = repmat(StimFrameIndex,numStims,1);
end

% Format index
if ~isempty(StimFrameIndex)
    StimFrameIndex = bsxfun(@plus,[0;cumsum(numFrames(1:end-1))],StimFrameIndex); % offset by cumulative stimuli before each stimulus
    temp = false(1,size(Images{1},N));
    for sindex = 1:numStims
        if sindex ~= ControlIndex
            temp(StimFrameIndex(sindex,1):StimFrameIndex(sindex,2)) = true;
        end
    end
    StimFrameIndex = temp;
end


%% Determine text to display
if isequal(Text,true)
    Text = 1:numStims;
end
if ~isempty(Text)
    TextIndex = repelem(1:numStims,numFrames); % display text for each frame in that stimulus
end


%% Save to video
[Filename, CLim, indMap, Dim, Map] = vid(Filename, Images,...
    'StimFrameIndex',   StimFrameIndex,...
    'Text',             Text,...
    'TextIndex',        TextIndex,...
    'FontSize',         FontSize,...
    'Color',            Color,...
    'ColorBar',         showColorBar,...
    'CbLabel',          CbLabel,...
    'CbFontSize',       CbFontSize,...
    'Crop',             Crop,...
    'borderLims',       borderLims,...
    'Afilt',            Afilt,...
    'Tfilt',            Tfilt,...
    'Maps',             Maps,...
    'MergeType',        MergeType,...
    'Map',              indMap,Dim,Map,...
    'CMap',             CMap,...
    'CLim',             CLim,...
    'frameRate',        frameRate,...
    'pixelSize',        pixelSize,...
    'outputSize',       outputSize);


