function patchHandles = overlayROIs(Data, varargin)
% Data is a cell array of strings, cell array of vertices, or list of
% centroids

axesHandle = [];
Maps = [];
ROIindex = [1 inf];
FileIndex = [];

% Display Properties
Radius = [];
Color = [1,0,0];
EdgeColor = []; % defaults to 'Color'
LineWidth = 1;
FaceAlpha = 1;
EdgeAlpha = 1;
FaceBrightness = 1;
EdgeBrightness = 1;

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Axes', 'axes'}
                axesHandle = varargin{index+1};
                index = index + 2;
            case 'Maps'
                Maps = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FileIndex'
                FileIndex = varargin{index+1};
                index = index + 2;
            case 'Radius'
                Radius = varargin{index+1};
                index = index + 2;
            case 'Color'
                Color = varargin{index+1};
                index = index + 2;
            case 'EdgeColor'
                EdgeColor = varargin{index+1};
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
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Data', 'var') || isempty(Data)
    [Data, p] = uigetfile({'*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Data)
        return
    end
    Data = fullfile(p, Data);
end
if ischar(Data)
    Data = {Data};
end

if ~isempty(Radius)
    loadType = 'circle';
else
    loadType = 'roi';
end


%% Load ROIs
if iscellstr(Data) || isstruct(Data) || (iscell(Data) && isstruct(Data{1}))
    switch loadType
        case 'roi'
            [Data, FileIndex] = gatherROIdata(Data,'vertices',':','none',ROIindex,FileIndex,'outputFormat','cell');
            Data = cellfun(@(x) reshape(x, numel(x)/2, 2), Data, 'UniformOutput', false); % Reshape vertices
        case 'circle'
            [Data, FileIndex] = gatherROIdata(Data,'centroid',':','none',ROIindex,FileIndex);
    end
end

if isnumeric(Data)
    switch loadType
        case 'roi'
            Data = {Data};
        case 'circle'
            temp = Data;
            Data = cell(size(temp,1),1);
            for rindex = 1:size(temp,1)
                Data{rindex} = temp(rindex,:);
            end
    end
end

numROIs = numel(Data);

if isempty(FileIndex)
    FileIndex = ones(numROIs,1);
end


%% Build up Map and translate ROIs onto entire image
if ~isempty(Maps) % multiple files input

    % Load Maps
    if iscellstr(Maps)
        Files = Maps;
        Maps = imref2d;
        for findex = 1:numel(Files)
            load(Files{findex}, 'Map', '-mat');
            Maps(findex) = Map;
        end
    end
    
    % Build map
    offsets = mapFoVs(Maps, 'Type', 'index');
    
    % Translate ROIs
    for rindex = 1:numROIs
        Data{rindex} = bsxfun(@plus, Data{rindex}, offsets(FileIndex(rindex), 1:2));
    end
    
end


%% Determine properties
if size(Color,1)==1
    Color = repmat(Color,numROIs,1);
end
if isempty(EdgeColor)
    EdgeColor = Color;
elseif size(EdgeColor,1)==1
    EdgeColor = repmat(EdgeColor,numROIs,1);
end
if numel(LineWidth)==1
    LineWidth = repmat(LineWidth,numROIs,1);
end
if numel(FaceAlpha)==1
    FaceAlpha = repmat(FaceAlpha,numROIs,1);
end
if numel(EdgeAlpha)==1
    EdgeAlpha = repmat(EdgeAlpha,numROIs,1);
end
if numel(FaceBrightness)==1
    FaceBrightness = repmat(FaceBrightness,numROIs,1);
end
if numel(EdgeBrightness)==1
    EdgeBrightness = repmat(EdgeBrightness,numROIs,1);
end

if strcmp(loadType,'circle')
    if isempty(Radius)
        Radius = 2*ones(numROIs,1);
    elseif numel(Radius)==1
        Radius = repmat(Radius,numROIs,1);
    end
    for rindex = 1:numROIs
        Data{rindex} = circle(Data{rindex},Radius(rindex));
    end
end


%% Plot ROIs

% Create figure
if isempty(axesHandle)
    figure;
    axesHandle = gca;
end

% Plot each ROI
patchHandles = nan(numROIs,1);
for rindex = 1:numROIs
    patchHandles(rindex) = patch(...
        Data{rindex}(:,1),...
        Data{rindex}(:,2),...
        Color(rindex,:)*FaceBrightness(rindex),...
        'Parent',               axesHandle,...
        'FaceAlpha',            FaceAlpha(rindex),...
        'EdgeAlpha',            EdgeAlpha(rindex),...
        'EdgeColor',            EdgeColor(rindex,:)*EdgeBrightness(rindex),...
        'LineWidth',            LineWidth(rindex),...
        'UserData',             rindex);
    
end

end

function vertices = circle(centroid,r)
th = (0:pi/50:2*pi)';
vertices = [r * cos(th) + centroid(1), r * sin(th) + centroid(2)];
end
