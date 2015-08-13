function [Index, ROIMasks] = autoMatchROIs(ROIMasks, Maps, Centroids, DataIndex)

saveOut = false;
saveFile = {''};

distanceThreshold = 10; % pixels
overlapThreshold = .9; % percentage
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
numFiles = numel(ROIMasks);

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
    DataIndex = 1:numFiles;
end


%% Load in ROIMasks
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
[~,~,numROIs] = cellfun(@size, ROIMasks);

% Compute centroids (if not given as input argument or loaded)
if ~exist('Centroids', 'var') || isempty(Centroids)
    Centroids = cell(numFiles, 1);
    for findex = 1:numFiles
        Centroids{findex} = zeros(numROIs(findex), 2);
        for rindex = 1:numROIs(findex)
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


%% Determine location of images
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

% Determine image offsets
offsets = zeros(numFiles, 2);
for findex = 1:numFiles
    offsets(findex, :) = [Maps{findex}.XWorldLimits(1) - XLim(1), Maps{findex}.YWorldLimits(1) - YLim(1)];
end

% Build Map
if addNewROIs
    Map = zeros(H, W, numFiles);
    for findex = 1:numFiles
        ylim = [ceil(Maps{findex}.YWorldLimits(1)),floor(Maps{findex}.YWorldLimits(2))] - floor(YLim(1));
        xlim = [ceil(Maps{findex}.XWorldLimits(1)),floor(Maps{findex}.XWorldLimits(2))] - floor(XLim(1));
        Map(ylim(1):ylim(2),xlim(1):xlim(2),findex) = 1;
    end
end
Map = reshape(Map, H*W, numFiles);


%% Translate ROIs
ActualMasks = cell(numFiles,1);
ROIindex = cell(numFiles,1);
for findex = 1:numFiles
    ROIindex{findex} = 1:numROIs(findex);
    Centroids{findex} = bsxfun(@plus, Centroids{findex}, offsets(findex,:));
    ActualMasks{findex} = false(H, W, numROIs(findex));
    for rindex = 1:numROIs(findex)
        ActualMasks{findex}(:,:,rindex) = imwarp(ROIMasks{findex}(:,:,rindex), Maps{findex}, affine2d(), 'OutputView', MainMap);
    end
end
ActualMasks = cellfun(@(x) reshape(x, size(x,1)*size(x,2), size(x,3)), ActualMasks, 'UniformOutput', false);


%% Determine distances between ROI centroids in the different datasets
combinations = combnk(1:numFiles, 2);
[~,I] = sort(combinations, 1);
combinations = combinations(I(:,1), :); % sort order
numCombinations = size(combinations, 1);
distances = cell(numCombinations, 1);
for cindex = 1:numCombinations
    distances{cindex} = pdist2(Centroids{combinations(cindex, 1)}, Centroids{combinations(cindex, 2)});
end


%% Match ROIs
Index = [];
numUniqueROIs = 0;
for findex = 1:numFiles
    numCurrentROIs = numel(ROIindex{findex});               % unmatched ROIs left in current dataset
    Index = cat(1, Index, nan(numCurrentROIs, numFiles));   % expand list of unique ROIs
    Index(numUniqueROIs+1:end, findex) = ROIindex{findex};  % record ROI identifiers for current file
    
    % Cycle through remaining ROIs matching each to the other files
    for rindex = ROIindex{findex}
        
        % Determine what files current ROI should be found in
        FileIndices = find(all(Map(ActualMasks{findex}(:,rindex),:), 1));
        FileIndices(FileIndices==findex) = [];
        
        % Match ROIs in each of the matched Files
        for mfindex = FileIndices
            
            % Determine closest existing ROIs in matched file
            cindex = ismember(combinations,[findex, mfindex],'rows');
            if any(cindex)
                currentDistances = distances{cindex}(rindex, :);
            else
                cindex = ismember(combinations,[mfindex, findex],'rows');
                currentDistances = distances{cindex}(:, rindex);
            end
            [~,distIndices] = sort(currentDistances);
            
            % Determine if any of the ROIs within the distance threshold
            % overlap more than the overlap threshold
            for mrindex = distIndices
                distThresh = currentDistances(mrindex) <= distanceThreshold;
                overlapThresh = nnz(ActualMasks{findex}(:,rindex) & ActualMasks{mfindex}(:,mrindex))/nnz(ActualMasks{findex}(:,rindex)) >= overlapThreshold;
                if ~distThresh                                          % moved onto ROIs too far away
                    Index(numUniqueROIs+rindex, mfindex) = 0;           % therefore ROI doesn't exist in matched file
                    break
                elseif overlapThresh
                    Index(numUniqueROIs+rindex, mfindex) = mrindex;     % ROI matches
                    ROIindex{mfindex}(ROIindex{mfindex}==mrindex) = []; % remove matching ROI
                    break
                end
            end %mrindex
        end %mfindex
        
        % Remove current ROI so can't be matched with other ROIs in future
        ROIindex{findex}(ROIindex{findex}==rindex) = []; % remove current ROI from list
        
    end %rindex
    
    numUniqueROIs = numUniqueROIs + numCurrentROIs; % update ROI counter
    
end % findex


%% Add unmatched ROIs
if addNewROIs
    [uindex, findex] = find(Index == 0);
    for rindex = 1:numel(uindex)
        [mfindex, ~, mrindex] = find(Index(uindex(rindex),:), 1);                                           % find matched ROI in another file
        newMask = imwarp(ROIMasks{mfindex}(:,:,mrindex), Map{mfindex}, 'OutputView', Map{findex(rindex)});  % translate ROI
        ROIMasks{findex(rindex)} = cat(3, ROIMasks{findex(rindex)}, newMask);                               % add new ROI
        Index(uindex(rindex), findex(rindex)) = size(ROIMasks{findex(rindex)}, 3);                          % record index of new ROI
    end
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