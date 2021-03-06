function [img, BW, Ind] = embedROIs(Vertices, varargin)

showImg = false;

Color = 1;
CMap = [];
hA = [];
Maps = [];
img = [];
FileIndex = [];

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Axes', 'axes'}
                hA = varargin{index+1};
                index = index + 2;
            case 'Maps'
                Maps = varargin{index+1};
                index = index + 2;
            case 'FileIndex'
                FileIndex = varargin{index+1};
                index = index + 2;
            case 'Color'
                Color = varargin{index+1};
                index = index + 2;
            case {'ColorMap','colormap','cmap','CMap','Colormap'}
                CMap = varargin{index+1};
                index = index + 2;
            case {'Image','image','img'}
                img = varargin{index+1};
                index = index + 2;
            case 'show'
                showImg = true;
                index = index + 1;
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
    [Vertices, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Vertices)
        return
    end
    Vertices = fullfile(p, Vertices);
end
if ~iscell(Vertices) % single ROI input
    Vertices = {Vertices};
end


%% Load ROIs
if iscellstr(Vertices) || isstruct(Vertices) || (iscell(Vertices) && isstruct(Vertices{1}))
    [Vertices, FileIndex] = gatherROIdata(Vertices, 'vertices', ':', 'none');
    Vertices = cellfun(@(x) reshape(x, numel(x)/2, 2), Vertices, 'UniformOutput', false); % Reshape vertices
end
numROIs = numel(Vertices);


%% mFoV: Build up Map and translate ROIs onto full image
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
    [offsets, refMap] = mapFoVs(Maps, 'Type', 'index');
    
    % Translate ROIs
    for rindex = 1:numROIs
        Vertices{rindex} = bsxfun(@plus, Vertices{rindex}, offsets(FileIndex(rindex),1:2));
    end
    
else
    refMap = [];    
end


%% Create image

% Determine image dimensions
if ~isempty(img)
    [H,W] = size(img);
else
    if ~isempty(refMap)
        H = refMap.ImageExtentInWorldY;
        W = refMap.ImageExtentInWorldX;
    else
        H = 512;
        W = 796;
    end
    img = nan(H,W); % create image
end

% Determine ROI color
if numel(Color)==1
    Color = repmat(Color, numROIs, 1);
end

% Display each ROI
BW = false(size(img,1),size(img,2),numROIs);
for rindex = numROIs:-1:1 % ensure early ROIs are on top
    BW(:,:,rindex) = poly2mask(Vertices{rindex}(:,1),Vertices{rindex}(:,2),H,W);
    img(BW(:,:,rindex)) = Color(rindex);
end

% Convert large binary matrix to indexed array
if nargout>2
    temp = cat(3,false(H,W),BW);
    [~,Ind] = max(temp,[],3);
    Ind = Ind - 1;
end


%% Display image
if showImg
    
    % Determine colormap
    if isempty(CMap)
        % CMap = b2r(min(img(:)), max(img(:)));
        CMap = HiLoColormap(min(img(:)), max(img(:)));
    end
    
    % Create figure
    if isempty(hA)
        figure;
        hA = gca;
    end
    axes(hA);
    
    % Display image
    image(img);
    colormap(CMap);
    
end
