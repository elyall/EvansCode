function [Centroids, Vertices, ROIMasks] = mapROIs(ROIs, Maps)



%% Parse input arguments
if ~exist('ROIs', 'var')
    [ROIs,p] = uigetfile({'*.rois;*.mat'},'Choose ROI files','MultiSelect','on');
    if isnumeric(ROIs)
        return
    elseif ischar(ROIs)
        ROIs = {fullfile(p,ROIs)};
    else
        ROIs = fullfile(p, ROIs);
    end
end

if ~exist('Maps', 'var')
    [Maps,p] = uigetfile({'*.exp;*.align'}, 'Select files containing maps:', directory, 'MultiSelect', 'on');
    if isnumeric(Maps)
        return
    elseif iscellstr(Maps)
        Maps = fullfile(p, Maps);
    else
        Maps = {fullfile(p, Maps)};
    end
end

numFiles = numel(ROIs);


%% Load in data

% Load in ROIs
if iscellstr(ROIs)
    ROIFiles = ROIs;
    ROIs = cell(numFiles, 1);
    for findex = 1:numFiles
        load(ROIFiles{findex}, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata;
    end
end
numROIs = cellfun(@(x) numel(x.rois), ROIs);

% Load in maps
if iscellstr(Maps)
    MapFiles = Maps;
    Maps = cell(numFiles, 1);
    for findex = 1:numFiles
        load(MapFiles{findex}, 'Map', '-mat');
        Maps{findex} = Map;
    end
end


%% Determine file offsets
[Offsets, Map, XLim, YLim] = mapFoVs(Maps); % determine amount of translation between the datasets


%% Translate ROI centroids
Centroids = cell(numFiles, 1); % initialize output
for findex = 1:numFiles
    Centroids{findex} = reshape([ROIs{findex}.rois(:).centroid], 2, numROIs(findex))';  % gather centroids
    Centroids{findex} = bsxfun(@plus, Centroids{findex}, Offsets(findex,:));            % translate centroids
end


%% Translate ROI vertices
if nargout > 1
    Vertices = cell(numFiles, 1); % initialize output
    for findex = 1:numFiles
        Vertices{findex} = cell(numROIs(findex), 1); % initialize for 'deal'
        [Vertices{findex}{:}] = deal(ROIs{findex}.rois.vertices); % gather vertices
        Vertices{findex} = cellfun(@(x) bsxfun(@plus, x, Offsets(findex,:)), Vertices{findex}, 'UniformOutput', false); % translate vertices
    end
end


%% Create mFoV ROI masks
if nargout > 2
    
    % Initialize reference
    [H,W,~] = size(Map);
    MainMap = imref2d([H,W], XLim, YLim);
    
    % Create new masks
    ROIMasks = cell(numFiles,1); % initialize output
    for findex = 1:numFiles
        ROIMasks{findex} = false(H, W, numROIs(findex));
        for rindex = 1:numROIs(findex)
            ROIMasks{findex}(:,:,rindex) = imwarp(ROIs{findex}.rois(rindex).mask, Maps{findex}, affine2d(), 'OutputView', MainMap);
        end
    end
end


