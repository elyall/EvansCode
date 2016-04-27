function [Centroids, Vertices, ROIMasks] = mapROIs(ROIs, Maps)

% ROIindex = [1 inf];
% FileIndex = [1 inf];


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
    Maps = imref2d();
    for findex = 1:numFiles
        temp = load(MapFiles{findex}, 'Map', '-mat');
        if isfield(temp,'Map')
            Maps(findex) = temp.Map;
        end
        clear temp;
    end
end


%% Determine file offsets
[Offsets, refMap] = mapFoVs(Maps); % determine amount of translation between the datasets


%% Translate ROI centroids
Centroids = cell(numFiles, 1); % initialize output
for findex = 1:numFiles
    Centroids{findex} = reshape([ROIs{findex}.rois(:).centroid], 2, numROIs(findex))';  % gather centroids
    Centroids{findex} = bsxfun(@plus, Centroids{findex}, Offsets(findex,[1,2]));            % translate centroids
end


%% Translate ROI vertices
if nargout > 1
    Vertices = cell(numFiles, 1); % initialize output
    for findex = 1:numFiles
        Vertices{findex} = cell(numROIs(findex), 1); % initialize for 'deal'
        [Vertices{findex}{:}] = deal(ROIs{findex}.rois.vertices); % gather vertices
        Vertices{findex} = cellfun(@(x) bsxfun(@plus, x, Offsets(findex,[1,2])), Vertices{findex}, 'UniformOutput', false); % translate vertices
    end
end


%% Create mFoV ROI masks
if nargout > 2
    ROIMasks = cell(numFiles,1); % initialize output
    for findex = 1:numFiles
        ROIMasks{findex} = false(refMap.ImageExtentInWorldY, refMap.ImageExtentInWorldX, numROIs(findex));
        for rindex = 1:numROIs(findex)
            ROIMasks{findex}(:,:,rindex) = imwarp(ROIs{findex}.rois(rindex).mask, Maps{findex}, affine2d(), 'OutputView', refMap);
        end
    end
end


