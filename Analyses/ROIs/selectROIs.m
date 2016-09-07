function [inside, ROIindex, FileIndex, vertices] = selectROIs(Data, vertices, varargin)

ROIindex = [];  % indices of ROIs input (or should load)
FileIndex = []; % indices of ROIs input (or should load)
Maps = [];      % for multiple files
Image = [];     % for UI selection

timeToEditROI = 5;  % seconds
verbose = false;    % display figure at end
directory = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FileIndex'
                FileIndex = varargin{index+1};
                index = index + 2;
            case 'Maps'
                Maps = varargin{index+1};
                index = index + 2;
            case 'Image'
                Image = varargin{index+1};
                index = index + 2;
            case {'time','timeToEditROI'}
                timeToEditROI = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
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

if isequal(Data,true)
    [Data, p] = uigetfile({'*.rois'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Data)
        return
    end
    Data = fullfile(p,Data);
end

if ~exist('vertices','var') || isempty(vertices)
    vertices = [];
end

if isequal(Maps,true)
    [Maps, p] = uigetfile({'*.mat;*.ciexp'},'Choose Experiment file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Maps)
        return
    end
    Maps = fullfile(p,Maps);
end


%% Load in data
if iscellstr(Data) || isstruct(Data) || (iscell(Data) && isstruct(Data{1}))
    [Data, FileIndex, ROIindex] = gatherROIdata(Data, 'centroid', ':', 'none', ROIindex, FileIndex);
else
    if isempty(FileIndex)
        FileIndex = (1:size(Data,1))';
    end
    if isempty(ROIindex)
        ROIindex = nan(numel(FileIndex),1);
        for findex = unique(FileIndex)'
            ROIindex(FileIndex==findex) = 1:nnz(FileIndex==findex);
        end
    end
end


%% Build up Map and translate ROIs onto entire image
if ~isempty(Maps) % multiple files input

    % Load Maps
    if ischar(Maps)
        Maps = {Maps};
    end
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
    if isnumeric(Data) % centroids
        for rindex = unique(FileIndex)'
            Data(FileIndex==rindex,:) = bsxfun(@plus, Data(FileIndex==rindex,:), offsets(rindex,1:2));
        end
    elseif iscell(Data) % vertices
        for rindex = 1:numel(Data)
            Data{rindex} = bsxfun(@plus, Data{rindex}, offsets(FileIndex(rindex),1:2));
        end
    end
end


%% Select ROIs
if isempty(vertices) || isequal(vertices, true)
    hF = figure;
    if ~isempty(Image)
        imagesc(Image); hold on;
    end
    if ~isempty(Data)
        if isnumeric(Data)
            plot(Data(:,1),Data(:,2),'x.');
        elseif iscell(Data)
            for rindex = 1:numel(Data)
                patch(Data{rindex}(:,1),Data{rindex}(:,2),'k');
            end
        end
    end
    if isempty(Image)
        set(gca,'YDir','reverse');
    end
    h = impoly(gca,'Closed',true);
    api = iptgetapi(h);
    fcn = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
    api.setPositionConstraintFcn(fcn);
    for index = timeToEditROI:-1:1
        set(hF, 'Name', sprintf('Closing in %d seconds', index));
        pause(1);
    end
    vertices = getPosition(h);
    vertices = [vertices; vertices(1,:)];
    close(gcf);
end


%% Determine which ROIs are inside the selected region
if ~isempty(Data)
    if isnumeric(Data)  % centroids
        inside = logical(inpolygon(Data(:,1),Data(:,2),vertices(:,1),vertices(:,2)));
    elseif iscell(Data) % vertices
        inside = false(numel(Data),1);
        for rindex = 1:numel(Data)
            inside(rindex) = all(inpolygon(Data{rindex}(:,1),Data{rindex}(:,2),vertices(:,1),vertices(:,2)));
        end
    end
    ROIindex = ROIindex(inside);
    FileIndex = FileIndex(inside);
else
    inside = [];
    ROIindex = [];
    FileIndex = [];
end


%% Display results
if verbose
    figure; hold on;
    if ~isempty(Image)
        imagesc(Image); hold on;
    end
    if ~isempty(Data)
        if isnumeric(Data)  % centroids
            plot(Data(inside,1),Data(inside,2),'rx');
            plot(Data(~inside,1),Data(~inside,2),'bx');
        elseif iscell(Data) % vertices
            for rindex = 1:numel(Data)
                if inside(rindex)
                    patch(Data{rindex}(:,1),Data{rindex}(:,2),'r');
                else
                    patch(Data{rindex}(:,1),Data{rindex}(:,2),'b');
                end
            end
        end
    end
    plot(vertices(:,1),vertices(:,2),'g-');
    if isempty(Image)
        set(gca,'YDir','reverse');
    end
end

