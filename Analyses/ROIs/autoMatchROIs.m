function [Index, ROIMasks] = autoMatchROIs(ROIMasks, Maps, Centroids, DataIndex)

saveOut = false;

distanceThreshold = 10; % pixels
overlapThreshold = 90; % percentage
addNewROIs = true; % add ROIs that don't match across files

directory = cd;

%% Parse input arguments
if ~exist('ROIMasks', 'var') || isempty(ROIMasks)
    [ROIMasks,p] = uigetfile({'*.rois;*.mat'}, 'Select ROI files:', directory, 'MultiSelect', 'on');
    if isnumeric(ROIMasks)
        return
    elseif iscellstr(ROIMasks)
        ROIMasks = fullfile(p, ROIMasks);
    else
        ROIMasks = {fullfile(p, ROIMasks)};
    end
end

if ~exist('Maps', 'var') || isempty(Maps)
    [Maps,p] = uigetfile({'*.exp;*.align'}, 'Select files containing maps:', directory, 'MultiSelect', 'on');
    if isnumeric(Maps)
        return
    elseif iscellstr(Maps)
        Maps = fullfile(p, Maps);
    else
        Maps = {fullfile(p, Maps)};
    end
end

if ~exist('DataIndex', 'var') || isempty(DataIndex)
    DataIndex = 1:numel(ROIMasks);
end


%% Load in ROIMasks
numFiles = numel(ROIMasks);
if iscellstr(ROIMasks)
    ROIFiles = ROIMasks;
    InitialROIdata = cell(numFiles, 1);
    ROIMasks = cell(numFiles, 1);
    Centroids = cell(numFiles, 1);
    for findex = 1:numFiles
        [~,~,ext] = fileparts(ROIFiles{findex});
        switch ext
            case '.segment'
                load(ROIFiles{findex}, 'mask', '-mat');
                ROIMasks{findex} = mask;
                Centroids{findex} = zeros(size(ROIMasks{findex},3), 2);
                for rindex = 1:size(ROIMasks{findex},3)
                    temp = regionprops(ROIMasks{findex}(:,:,rindex), 'centroid');
                    Centroids{findex}(rindex,:) = temp.Centroid;
                end
            case '.rois'
                load(ROIFiles{findex}, 'ROIdata', '-mat');
                InitialROIdata{findex} = ROIdata;
                ROIMasks{findex} = reshape([ROIdata.rois(:).pixels], size(ROIdata.rois(1).pixels,1), size(ROIdata.rois(1).pixels,2), numel(ROIdata.rois));
                Centroids{findex} = reshape([ROIdata.rois(:).centroid], 2, numel(ROIdata.rois))';
        end
    end
end

% Compute centroids (if not given as input argument or loaded)
if ~exist('Centroids', 'var') || isempty(Centroids)
    Centroids = cell(numFiles, 1);
    for findex = 1:numFiles
        Centroids{findex} = zeros(size(ROIMasks{findex},3), 2);
        for rindex = 1:size(ROIMasks{findex},3)
            temp = regionprops(ROIMasks{findex}(:,:,rindex), 'centroid');
            Centroids{findex}(rindex,:) = temp.Centroid;
        end
    end
end


%% Load in Maps
if iscellstr(Maps)
    MapFiles = Maps;
    Maps = cell(numel(MapFiles), 1);
    for findex = 1:numel(MapFiles)
        load(MapFiles{findex}, 'Map', '-mat');
        Maps{findex} = Map;
    end
end


%% Determine composite image size
XLim = [inf, -inf];
YLim = [inf, -inf];
for findex = 1:numFiles
    XLim(1) = min(XLim(1), Maps{findex}.XWorldLimits(1));
    XLim(2) = max(XLim(2), Maps{findex}.XWorldLimits(2));
    YLim(1) = min(YLim(1), Maps{findex}.YWorldLimits(1));
    YLim(2) = max(YLim(2), Maps{findex}.YWorldLimits(2));
end
H = diff(YLim);
W = diff(XLim);
MainMap = imref2d([H,W], XLim, YLim);


%% Determine offsets
offsets = zeros(numFiles, 2);
for findex = 1:numFiles
    offsets = [Maps{findex}.XWorldLimits(1) - XLim(1), Maps{findex}.YWorldLimits(1) - YLim(1)];
end


%% Translate ROIs
ActualMasks = cell(numFiles,1);
for findex = 1:numFiles
    ActualMasks{findex} = zeros(H, W, size(ROIMasks{findex}, 3));
    for rindex = 1:size(ROIMasks, 3)
        Centroids{findex}(rindex, :) = Centroids{findex}(rindex, :) + offsets(findex,:);
        ActualMasks{findex}(:,:,rindex) = imwarp(ROIMasks{findex}(:,:,rindex), Maps{findex}, affine2d(), 'OutputView', MainMap);
    end
end


%% Compute distance between all ROIs
distances = squareform(pdist(centers));


%% Match ROIs that fall within threshold from each other
% BaseTag = strcat(ID,'_',datestr(now,'yyyy-mm-dd'),'_');
% TagIndex = 1;
for rindex = 1:numROIs-1
    
    % Find close enough ROIs
    MatchedROIs = find(distances(rindex,rindex+1:end) <= distanceThreshold)+rindex;
    
    % Remove multiple ROIs from original file
    FileIndex = ROIs.rois(rindex).index(1);
    MatchedInfo = cat(2, distances(rindex,MatchedROIs)', reshape([ROIs.rois(MatchedROIs).index], 2, numel(MatchedROIs))');
    bad = MatchedInfo(:,2) == FileIndex;
    MatchedROIs(bad) = [];
    MatchedInfo(bad,:) = [];
    
    % Remove multiple ROIs from any file
    [~,index] = sort(MatchedInfo(:,1));
    MatchedROIs = MatchedROIs(index);
    MatchedInfo = MatchedInfo(index, :);
    for findex = 2:numFiles
        index = MatchedInfo(:,2) == findex;
        nMatchesFromFile = sum(index);
        if nMatchesFromFile > 1
            bad = find(index, nMatchesFromFile-1, 'last'); % keep only closest match
            MatchedROIs(bad) = [];
            MatchedInfo(bad, :) = [];
        end
    end
    numrois = numel(MatchedROIs);
    
%     % Determine if any matched ROIs have already been assigned a tag
%     hasTag = false(numrois+1, 1);
%     if ~all(ismember(ROIs.rois(rindex).tag,'0123456789')); % already assigned tag
%         hasTag = true;
%     end
%     for mindex = 1:numrois
%         if ~all(ismember(ROIs.rois(MatchedROIs(mindex)).tag,'0123456789'));
%                 hasTag(mindex+1) = true;
%         end
%     end
    
%     % Determine or assign tag
%     if sum(hasTag) > 1
%         warning('Overwriting previously assigned tag');
%     end
%     if any(hasTag) % previously assigned tag
%         mindex = find(hasTag, 1);
%         if mindex == 1
%             Tag = ROIs.rois(rindex).tag;
%         else
%             Tag = ROIs.rois(MatchedROIs(mindex-1)).tag;
%         end
%     else % assign new tag
%         Tag = strcat(BaseTag, num2str(TagIndex));
%         TagIndex = TagIndex + 1;
%     end
 
%     % Save Tags
%     Data(ROIs.rois(rindex).index(1)).ROIdata.rois(ROIs.rois(rindex).index(2)).tag = Tag;
%     ROIs.rois(rindex).Tag = Tag;
%     for mindex = 1:numrois
%         Data(ROIs.rois(MatchedROIs(mindex)).index(1)).ROIdata.rois(ROIs.rois(MatchedROIs(mindex)).index(2)).tag = Tag;
%         ROIs.rois(MatchedROIs(mindex)).tag = Tag;
%     end
end

%% Save output
if saveOut && ~isempty(saveFile)
    for findex = 1:numel(saveFile)
        [~,~,ext] = fileparts(saveFile{findex});
        switch ext
            case '.segment'
                mask = ROIMasks{findex};
                if ~exist(saveFile{findex}, 'file')
                    save(saveFile{findex}, 'mask', '-mat', '-v7.3');
                else
                    save(saveFile{findex}, 'mask', '-mat', '-append');
                end
            case '.rois'
                ROIdata = createROIdata(ROIMasks{findex}(:,:,newMasks:end), 'ROIdata', InitialROIdata{findex});
                if ~exist(saveFile{findex}, 'file')
                    save(saveFile{findex}, 'ROIdata', '-mat', '-v7.3');
                else
                    save(saveFile{findex}, 'ROIdata', '-mat', '-append');
                end
        end
    end
end