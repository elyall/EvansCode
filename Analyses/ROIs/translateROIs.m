function ROIs = translateROIs(ROIs, XYShifts, Maps)

%% Parse input argumetns
if ~exist('ROIs','var') || isempty(ROIs)
    [ROIs,p] = uigetfile({'*.rois;*.mat'}, 'Select ROI files:', directory, 'MultiSelect', 'on');
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p,ROIs);
end
if ischar(ROIs) || isstruct(ROIs)
    ROIs = {ROIs};
end
numFiles = numel(ROIs);

if isa(XYShifts, 'imref2d')
    temp = zeros(numel(XYShifts),2);
    for findex = 1:numel(XYShifts)
        temp(findex,:) = [XYShifts.YWorldLimits(1), XYShifts.XWorldLimits(1)] - .5;
    end
    XYShifts = temp;
end
if size(XYShifts,1)~=numFiles
    XYShifts = repmat(XYShifts,numFiles,1);
end


%% Load in data
if iscellstr(ROIs)
    for findex = 1:numFiles
        temp = load(ROIs{findex},'ROIdata','-mat');
        ROIs{findex} = temp.ROIdata;
    end
end

if ~exist('Maps','var') || isempty(Maps)
    Maps = imref2d();
    for findex = 1:numFiles
        Maps(findex) = imref2d(size(ROIs{findex}.rois(1).pixels));
    end
end
if numel(Maps)~=numFiles
    Maps = repmat(Maps,1,numFiles);
end


%% Translate ROIs

% Translate pixels
for findex = 1:numFiles
    for rindex = 1:numel(ROIs{findex}.rois)
        ROIs{findex}.rois(rindex).vertices = bsxfun(@plus, ROIs{findex}.rois(rindex).vertices, XYShifts(findex,:));
    end
end

% Translate masks
YXShifts = flip(XYShifts,2);
for findex = 1:numFiles
    for rindex = 1:numel(ROIs{findex}.rois)
        ROIs{findex}.rois(rindex).pixels = sparse(circshift(full(ROIs{findex}.rois(rindex).pixels), YXShifts(findex,:)));
        ROIs{findex}.rois(rindex).mask = sparse(circshift(full(ROIs{findex}.rois(rindex).mask), YXShifts(findex,:)));
        ROIs{findex}.rois(rindex).neuropilmask = sparse(circshift(full(ROIs{findex}.rois(rindex).neuropilmask), YXShifts(findex,:)));
    end
end


