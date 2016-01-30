function patchHandles = overlayROIs(Vertices, varargin)

axesHandle = [];
Maps = [];
Color = 'm';
ObjSize = 1;

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
            case 'Color'
                Color = varargin{index+1};
                index = index + 2;
            case {'LineWidth','MarkerSize'}
                ObjSize = varargin{index+1};
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

if ~exist('Vertices', 'var') || isempty(Vertices)
    [Vertices, p] = uigetfile({'*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Vertices)
        return
    end
    Vertices = fullfile(p, Vertices);
end
if ~iscell(Vertices) % single ROI input
    Vertices = {Vertices};
end


%% Load ROIs
if iscellstr(Vertices)
    [Vertices, FileIndex] = gatherROIdata(Vertices, 'vertices', ':', 'none');
    Vertices = cellfun(@(x) reshape(x, numel(x)/2, 2), Vertices, 'UniformOutput', false); % Reshape vertices
end
numROIs = numel(Vertices);


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
        Vertices{rindex} = bsxfun(@plus, Vertices{rindex}, offsets(FileIndex(rindex), 1:2));
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
        Vertices{rindex}(:,1),...
        Vertices{rindex}(:,2),...
        Color,...
        'Parent',               axesHandle,...
        'FaceAlpha',            0,...
        'EdgeColor',            Color,...
        'LineWidth',            ObjSize,...
        'UserData',             rindex);
end

