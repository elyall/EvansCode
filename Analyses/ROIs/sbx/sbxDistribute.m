function ROIdata = sbxDistribute(fname, varargin)

saveOut = false;
ROIFile = '';

%% Parse input arguments
if ~exist('fname', 'var') || isempty(fname)
    [fname,p] = uigetfile({'*.segment'},'Select segment file:',directory);
    if isnumeric(fname)
        return
    end
    fname = fullfile(p, fname);
    fname = fname(1:end-4);
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
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


%% Initialize ROIdata structure
if isempty(ROIFile)
    ROIFile = [fname, '.rois'];
    if exist(ROIFile, 'file')
        index = 1;
        while exist([fname, '_', num2str(index), '.rois'], 'file')
            index = index + 1;
        end
        ROIFile = [fname, '_', num2str(index), '.rois'];
    end
end
ROIdata.filename = ROIFile;
ROIdata.offset = zeros(1,3);
ROIdata.imagefiles = {[fname,'.sbx']};
ROIdata.files = {[fname, '.sbx']};


%% Distribute ROI segmentation
load([fname, '.segment'], 'mask', '-mat');
numROIs = max(mask(:));
for rindex = 1:numROIs
    ROIdata.rois(rindex).mask = mask==rindex;
end


%% Filler information
for rindex = 1:numROIs
    ROIdata.rois(rindex).label = {};
    ROIdata.rois(rindex).tag = num2str(rindex);
    ROIdata.rois(rindex).color = [];
    ROIdata.rois(rindex).depth = 1;
    ROIdata.rois(rindex).frame = 1;
    ROIdata.rois(rindex).neuropilmask = [];
    
    ROIdata.rois(rindex).type = 'polygon';
    temp = bwboundaries(ROIdata.rois(rindex).mask);
    ROIdata.rois(rindex).vertices = flip(temp{1},2);
    ROIdata.rois(rindex).position = [];
    
    ROIdata.rois(rindex).pixels = ROIdata.rois(rindex).mask;
    
    temp = regionprops(ROIdata.rois(rindex).mask, 'Centroid');
    ROIdata.rois(rindex).centroid = temp.Centroid;
    
end


%% Distribute ROI signals
if exist([fname, '.signals'], 'file')
    load([fname, '.signals'], 'sig', 'pil', '-mat');
    for rindex = 1:numROIs
        ROIdata.rois(rindex).rawdata = sig(:,rindex)';
        ROIdata.rois(rindex).rawneuropil = pil(:,rindex)';
    end
end


%% Save ROIdata to file
if saveOut
    save(ROIFile, 'ROIdata', '-mat', '-v7.3');
end