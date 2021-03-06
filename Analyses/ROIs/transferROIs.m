function ROIMasks = transferROIs(ROIMasks, MapOrig, MapEnd, varargin)

saveOut = false;
saveFile = '';

ROIindex = [1 inf];

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{index+1};
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

if ~exist('ROIMasks', 'var') || isempty(ROIMasks)
    [ROIMasks,p] = uigetfile({'*.rois;*.segment'}, 'Select ROI file:', directory);
    if isnumeric(ROIMasks)
        return
    end
    ROIMasks = fullfile(p, ROIMasks);
end

if ~exist('MapOrig', 'var') || isempty(MapOrig)
    [MapOrig,p] = uigetfile({'*.exp;*.align'}, 'Select original Map file:', directory);
    if isnumeric(MapOrig)
        return
    end
    MapOrig = fullfile(p, MapOrig);
end

if ~exist('MapEnd', 'var') || isempty(MapEnd)
    [MapEnd,p] = uigetfile({'*.exp;*.align'}, 'Select end Map file:', directory);
    if isnumeric(MapEnd)
        return
    end
    MapEnd = fullfile(p, MapEnd);
end


%% Load in ROIs
if ischar(ROIMasks)
    ROIFile = ROIMasks;
    [~,~,ext] = fileparts(ROIFile);
    switch ext
        case '.rois'
            load(ROIFile, 'ROIdata', '-mat');
            ROIMasks = reshape(full([ROIdata.rois(:).pixels]), size(ROIdata.rois(1).pixels,1), size(ROIdata.rois(1).pixels,2), numel(ROIdata.rois));
        case '.segment'
            load(ROIFile, 'mask', 'dim', '-mat');
            if issparse(mask)
                ROIMasks = reshape(full(mask), dim(1), dim(2), size(mask,2));
            else
                ROIMasks = mask;
            end
    end
end

if numel(ROIindex)>1 && ROIindex(end)==inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(end-1)+1:size(ROIMasks,3));
end
ROIMasks = ROIMasks(:,:,ROIindex);
[H,W,numROIs] = size(ROIMasks);

% Report status
fprintf('Transferring %d ROI(s) ',numROIs);
if exist('ROIFile','var')
    fprintf('from %s ',ROIFile);
end
fprintf('to new axes...\t');


%% Load in Maps
if ischar(MapOrig)
    load(MapOrig, 'Map', '-mat');
    MapOrig = Map;
end
if ischar(MapEnd)
    if saveOut && isempty(saveFile) && exist('ext','var')
        [p,f,~] = fileparts(MapEnd);
        saveFile = fullfile(p,[f,ext]);
    end
    load(MapEnd, 'Map', '-mat');
    MapEnd = Map;
end


%% Create transformation object
x = MapOrig.XWorldLimits(1) - MapEnd.XWorldLimits(1);
y = MapOrig.YWorldLimits(1) - MapEnd.YWorldLimits(1);
tform = affine2d([1,0,0;0,1,0;x,y,1]);


%% Translate ROIs
for rindex = 1:numROIs
    ROIMasks(:,:,rindex) = imwarp(ROIMasks(:,:,rindex), MapOrig, tform, 'OutputView', MapOrig);
end
fprintf('Complete\n');


%% Remove ROIs that don't exist on new axes
emptyIndex = sum(reshape(ROIMasks, H*W, numROIs))==0;
if any(emptyIndex)
    fprintf('\t%d ROI(s) could not be transferred due to no overlap of FoVs: ',nnz(emptyIndex));
    fprintf('%d, ',find(emptyIndex));
    fprintf('\b\b\n');
    ROIMasks(:,:,emptyIndex) = [];
    numROIs = numROIs - nnz(emptyIndex);
end


%% Save ROIs to file
if saveOut && ~isempty(saveFile)
    [~,~,ext] = fileparts(saveFile);
    switch ext
        case '.segment'
            mask = sparse(reshape(ROIMasks, H*W, numROIs));
            if ~exist(saveFile, 'file')
                save(saveFile, 'mask', 'dim', '-mat', '-v7.3');
            else
                save(saveFile, 'mask', 'dim', '-mat', '-append');
            end
            fprintf('\tROI mask(s) saved to: %s\n',saveFile);
        case '.rois'
            ROIdata = createROIdata(ROIMasks, 'ROIdata', saveFile);     
            if ~exist(saveFile, 'file')
                save(saveFile, 'ROIdata', '-mat', '-v7.3');
            else
                save(saveFile, 'ROIdata', '-mat', '-append');
            end
            fprintf('\tROIdata generated and saved to: %s\n',saveFile);
    end
end

