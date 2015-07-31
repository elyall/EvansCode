function [ROIindex, FileIndex, Tags, ROIs] = determineROIs(ROIs, type)
% locate matched pre-trimming and post-trimming ROI

if ~exist('type', 'var')
    type = 'struct-index'; % 'tag' or 'struct-index'
end

%% Check input arguments
if ~exist('ROIs', 'var') || isempty(ROIs)
    directory = cd;
    [ROIs, p] = uigetfile({'*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(ROIs)
        return
    elseif iscell(ROIs)
        for findex = 1:numel(ROIs)
            ROIs{findex} = fullfile(p, ROIs{findex});
        end
    elseif ischar(ROIs)
        ROIs = {fullfile(p, ROIs)};
    end
end


%% Load data
numFiles = numel(ROIs);
if iscellstr(ROIs)
    ROIFiles = ROIs;
    ROIs = cell(numFiles, 1);
    for findex = 1:numFiles
        load(ROIFiles{findex}, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata;
    end
end

%% Pull Out Labels
switch type
    case 'tag' % match ROIs with the same Tag
        [~, fileindex, roiindex, ~, Tag, ~] = gatherROIdata(ROIs, 'rawdata', 1);
        [Tags, ~, ROIdict] = unique(Tag, 'stable');
        numROIs = numel(Tags);
        ROIindex = nan(numROIs, 2);
        FileIndex = nan(numROIs, 2);
        for rindex = 1:numROIs
            indices = find(ROIdict==rindex, 2);
            ROIindex(rindex, 1:numel(indices)) = roiindex(indices);
            FileIndex(rindex, 1:numel(indices)) = fileindex(indices);
        end
        
    case 'struct-index' % match ROIs by their index in the struct
        ROIindex = [];
        FileIndex = [];
        Tags = [];
        for findex = 1:2:numFiles
            numROIs = min(numel(ROIs{findex}.rois),numel(ROIs{findex+1}.rois));
            ROIindex = cat(1, ROIindex, repmat((1:numROIs)',1,2));
            FileIndex = cat(1, FileIndex, repmat([findex, findex+1],numROIs,1));
            Tags = cat(1, Tags, {ROIs{findex}.rois(1:numROIs).tag}');
        end
end
