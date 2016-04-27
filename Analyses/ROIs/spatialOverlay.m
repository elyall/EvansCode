function [Image, hA, patchHandles] = spatialOverlay(ROIs, Maps, ROIindex, FileIndex, ColorIndex, Labels, Colors, varargin)
%spatialOverlay    Overlay ROI data onto an image
% [Image, Origin, hA] = spatialOverlay(ROIs, Data, ROIindex, FileIndex, ColorIndex, Labels, Colors, Brightness)
%
% Labels - for discrete data, labels is a cell array of N strings
% corresponding to the <= N unique values in ColorIndex, and the N unique
% colors.
% Labels - for continuous data, labels is a 2x1 vector specifying the lower
% and upper bounds of the colors (the same thing as CLim)

% Image Properties
mergetype = 'quick'; % 'quick' or 'pretty'
Crop = false;
showColorBar = false;
colorbarLabel = '';
Title = '';

% ROI Properties
DataType = 'continuous'; % 'discrete' or 'continuous'
Radius = [];
LineWidth = 1;
FaceAlpha = 1;
EdgeAlpha = 1;
FaceBrightness = 1;
EdgeBrightness = 1;

% Default variables
Image = ones(512,796);
hA = [];
CMap = [];

% Miscellaneous (can't be changed via argument)
defaultColorMap = 'parula'; % 'parula' or 'jet'
FontSize_colorbarLabel = 20;
FontSize_colorbarTicks = 20;
flipColorbar = true;

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            
            % Image Properties
            case 'mergetype'
                mergetype = varargin{index+1};
                index = index + 2;
            case 'Crop'
                Crop = varargin{index+1};
                index = index + 2;
            case 'showColorBar'
                showColorBar = true;
                index = index + 1;
            case 'colorbarLabel'
                colorbarLabel = varargin{index+1};
                index = index + 2;
            case 'Title'
                Title = varargin{index+1};
                index = index + 2;
            
            % ROI Properties
            case 'DataType'
                DataType = varargin{index+1};
                index = index + 2;
            case 'Radius'
                Radius = varargin{index+1};
                index = index + 2;
            case 'LineWidth'
                LineWidth = varargin{index+1};
                index = index + 2;
            case 'FaceAlpha'
                FaceAlpha = varargin{index+1};
                index = index + 2;
            case 'EdgeAlpha'
                EdgeAlpha = varargin{index+1};
                index = index + 2;
            case 'FaceBrightness'
                FaceBrightness = varargin{index+1};
                index = index + 2;
            case 'EdgeBrightness'
                EdgeBrightness = varargin{index+1};
                index = index + 2;
            
            % Default variables
            case 'Image'
                Image = varargin{index+1};
                index = index + 2;
            case 'axes'
                hA = varargin{index+1};
                index = index + 2;
            case 'Colormap'
                CMap = varargin{index+1};
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

if ~exist('ROIs','var') || isempty(ROIs)
    [ROIs, p] = uigetfile({'*.rois'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p,ROIs);
end
if ischar(ROIs)
    ROIs = {ROIs};
end

if ~exist('Maps','var')
    [Maps, p] = uigetfile({'*.mat;*.ciexp'},'Choose Experiment file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Maps)
        return
    end
    Maps = fullfile(p,Maps);
end
if ischar(Maps)
    Maps = {Maps};
end


%% Load data

% Load ROIs
numFiles = numel(ROIs);
if iscellstr(ROIs)
    ROIFiles = ROIs;
    ROIs = cell(numFiles, 1);
    for findex = 1:numFiles
        load(ROIFiles{findex}, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata;
    end
end

% Load Maps
if iscellstr(Maps)
    MapFiles = Maps;
    Maps = imref2d();
    for findex = 1:numFiles
        temp = load(MapFiles{findex}, 'Map', '-mat');
        if isfield(temp,'Map')
            Maps(findex) = temp.Map;
        end
        clear temp;
    end
end


%% Display image

% Select axes
if isempty(hA)
    figure();
    hA = axes();
else
    axes(hA);
end

% Display Image
if ~isempty(Image)
        
    % Build Image
    if iscell(Image)
        if isempty(Maps); error('Maps required for creating image'); end
        Image = createImage(Image, Maps, 'Crop', Crop, 'speed', mergetype);
    end
    
    % Display Image
    if ~isempty(CMap)
        image(Image);
        colormap(CMap);
    else
        imagesc(Image)
    end
    axis equal off
    
end

% Display title
if ~isempty(Title)
    title(Title);
end


%% Determine ROIs to overlay
if ~exist('ROIindex', 'var') || isempty(ROIindex) || (ischar(ROIindex) && strcmp(ROIindex, 'all'))
    ROIindex = [1, inf];
end
if ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:numel(ROIs{1}.rois)); % defaults to only first file
end
numROIs = numel(ROIindex);

if ~exist('FileIndex', 'var') || isempty(FileIndex)
    FileIndex = ones(numROIs, 1); % defaults to only first file
end


%% Build ROI colormap
if ~isequal(round(ColorIndex),ColorIndex)
    DataType = 'continuous';
end

switch DataType
    case 'discrete'
        
        % Determine colors
        if ~exist('Colors', 'var') || isempty(Colors)
            N = numel(unique(ColorIndex));
            switch defaultColorMap
                case 'parula'
                    Colors = parula(N);
                case 'jet'
                    Colors = jet(N);
            end
        end
        
        % Determine labels
        if ~exist('Labels', 'var') || isempty(Labels)
            Labels = cellstr(num2str(unique(ColorIndex)));
            % Labels = cellstr(num2str((1:size(Colors,1))'));
        end
        
    case 'continuous'
        
        % Determine labels
        if ~exist('Labels', 'var') || isempty(Labels)
            Labels = [min(ColorIndex), max(ColorIndex)];
        end
        if numel(Labels)==2
            Labels = min(Labels):range(Labels)/9:max(Labels);
        end
        if isnumeric(Labels)
            if isrow(Labels)
                Labels = Labels';
            end
            Labels = cellstr(num2str(Labels));
        end
        
        % Change continuous space to discrete space
        if ~isequal(round(ColorIndex),ColorIndex)
            [ColorIndex, Labels, numSamples] = scaleContinuousData(ColorIndex);
        else
            numSamples = 255;
        end
        ColorIndex = ColorIndex - min(ColorIndex) + 1; % shift bottom to 1
        
        % Determine colors
        switch defaultColorMap
            case 'parula'
                Colors = parula(numSamples);
            case 'jet'
                Colors = jet(numSamples);
        end
        
end

if flipColorbar
    Labels = flip(Labels);
    Colors = flip(Colors);
    ColorIndex = abs(ColorIndex - max(ColorIndex) - 1);
end


%% Overlay ROIs
patchHandles = overlayROIs(ROIs,...
    'axes', hA, 'Maps', Maps, 'ROIindex', ROIindex, 'FileIndex', FileIndex,...
    'Radius', Radius, 'Color', Colors(ColorIndex,:), 'LineWidth', LineWidth,...
    'FaceAlpha', FaceAlpha, 'EdgeAlpha', EdgeAlpha, 'FaceBrightness', FaceBrightness,...
    'EdgeBrightness', EdgeBrightness);


%% Plot colorbar
if showColorBar
    
    % Determine colorbar properties
    cmap = colormap;                % determine colormap of image
    cbH = colorbar;                 % place colorbar
    YLim = get(cbH, 'Limits');      % determine colorbar limits
    NewCMap = cat(1,cmap,Colors);   % concatenate plot colors to colormap
    colormap(NewCMap);              % set new colormap
    HeightNewCMapInColorbar = size(Colors,1)*(YLim(2)/size(NewCMap,1)); % determine portion of colorbar taken by colors that were added
    NewYLim = [YLim(2)-HeightNewCMapInColorbar, YLim(2)];               % limit the colormap to this range
    Split = range(NewYLim)/numel(Labels);                               % determine the distance on the colorbar between two colors
    
    % Display labels
    switch DataType
        case 'discrete'
            set(cbH, 'Limits', NewYLim, 'FontSize', FontSize_colorbarTicks, 'Ticks', NewYLim(1)+Split/2:Split:NewYLim(2), 'YTickLabel', Labels);
        case 'continuous'
            set(cbH, 'Limits', NewYLim, 'FontSize', FontSize_colorbarTicks, 'Ticks', NewYLim(1):range(NewYLim)/(numel(Labels)-1):NewYLim(2), 'YTickLabel', Labels);
    end
    
    % Display colorbar title
    if ~isempty(colorbarLabel)
        ylabel(cbH, colorbarLabel, 'FontSize', FontSize_colorbarLabel);
    end
end


