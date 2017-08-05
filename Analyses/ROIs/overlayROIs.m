function [Handles,hA] = overlayROIs(Data, varargin)
% Data is a cell array of strings, cell array of vertices, or list of
% centroids

hA = [];
Maps = [];
ROIindex = [1 inf];
FileIndex = [];

% Display Properties
roiType = '';           % 'roi', 'circle', or 'arrow'
Radius = [];            % ('circle' or 'arrow')
FaceColor = [1,0,0];
EdgeColor = [];         % defaults to 'Color'
LineWidth = 1;
FaceAlpha = 1;
EdgeAlpha = 1;
FaceBrightness = 1;
EdgeBrightness = 1;
theta = []; % arrow only

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Axes','axes','axis','hA'}
                hA = varargin{index+1};
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
            case {'roiType','type','Type'}
                roiType = varargin{index+1};
                index = index + 2;
            case 'Radius'
                Radius = varargin{index+1};
                index = index + 2;
            case {'Color','FaceColor'}
                FaceColor = varargin{index+1};
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
            case 'theta'
                theta = varargin{index+1};
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
    [Data, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Data)
        return
    end
    Data = fullfile(p, Data);
end
if ischar(Data)
    Data = {Data};
end

if isempty(roiType) && ~isempty(Radius)
    roiType = 'circle';
elseif isempty(roiType) && isempty(Radius)
    roiType = 'roi';
end


%% Load ROIs
if iscellstr(Data) || isstruct(Data) || (iscell(Data) && isstruct(Data{1}))
    switch roiType
        case 'roi'
            [Data, FileIndex] = gatherROIdata(Data,'vertices',':','none',ROIindex,FileIndex,'outputFormat','cell');
            Data = cellfun(@(x) reshape(x, numel(x)/2, 2), Data, 'UniformOutput', false); % Reshape vertices
            ROIindex = 1:numel(Data);
        case {'circle','arrow'}
            [Data, FileIndex] = gatherROIdata(Data,'centroid',':','none',ROIindex,FileIndex);
            ROIindex = 1:size(Data,1);
    end
end

if isnumeric(Data)
    switch roiType
        case 'roi' % assumes one ROI
            Data = {Data};
        case {'circle','arrow'} % convert matrix to cell array
            Data = mat2cell(Data,ones(size(Data,1),1),2);
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
if size(FaceColor,1)==1
    FaceColor = repmat(FaceColor,numROIs,1);
end
if isempty(EdgeColor)
    EdgeColor = FaceColor;
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

if ismember(roiType,{'circle','arrow'}) && isempty(Radius)
    Radius = 2*ones(numROIs,1);
elseif ismember(roiType,{'circle','arrow'}) && numel(Radius)==1
    Radius = repmat(Radius,numROIs,1);
end

if strcmp(roiType,'circle') % create circle from centroids
    for rindex = 1:numROIs
        Data{rindex} = circle(Data{rindex},Radius(rindex)); %subfunction below
    end
end

if strcmp(roiType,'roi') % append starting point to end to complete polygon
    Data = cellfun(@(x) x([1:end,1],:),Data,'UniformOutput',false);
end

if strcmp(roiType,'arrow') % create arrows
    [x,y] = pol2cart(theta*pi/180,1);
    Radius = bsxfun(@times, repmat([x,y],numROIs,1), Radius);
    if iscell(Data)
        Data = cat(1,Data{:});
    end
end


%% Plot ROIs

% Create figure
if isempty(hA)
    figure;
    hA = gca;
end

% Plot each ROI
switch roiType
    
    case {'roi','circle'}
        Handles = nan(numROIs,1);
        for rindex = ROIindex
            Handles(rindex) = patch(...
                Data{rindex}(:,1),...
                Data{rindex}(:,2),...
                FaceColor(rindex,:)*FaceBrightness(rindex),...
                'Parent',               hA,...
                'FaceAlpha',            FaceAlpha(rindex),...
                'EdgeAlpha',            EdgeAlpha(rindex),...
                'EdgeColor',            EdgeColor(rindex,:)*EdgeBrightness(rindex),...
                'LineWidth',            LineWidth(rindex),...
                'UserData',             rindex);
        end
        
    case 'arrow'
        Handles = nan(numROIs,1);
        
        % Script from mathworks forums
        for rindex = ROIindex
            Handles(rindex) = arrow(Data(rindex,:),Data(rindex,:)+Radius(rindex,:),...
                'FaceColor',    FaceColor(rindex,:)*FaceBrightness(rindex),...
                'EdgeColor',    EdgeColor(rindex,:)*EdgeBrightness(rindex),...
                'Width',        LineWidth(rindex),...
                'BaseAngle',    60,...
                'TipAngle',     20,...
                'Length',       6);
        end
        
        % Quiver method
%         Handles = quiver(Data(:,1),Data(:,2),Radius(:,1),Radius(:,2),...
%             'Color',                    FaceColor(1,:)*FaceBrightness(1),...
%             'LineWidth',                LineWidth(1));
end

end

function vertices = circle(centroid,r)
th = (0:pi/50:2*pi)';
vertices = [r * cos(th) + centroid(1), r * sin(th) + centroid(2)];
end
