function ROIdata = computeDFoF(ROIdata, NeuropilWeight, varargin)

saveOut = false;
saveFile = '';
numFramesBaseline = 20;

%% Check input arguments
if ~exist('ROIdata','var') || isempty(ROIdata)
    directory = CanalSettings('DataDirectory');
    [ROIdata, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if ~exist('NeuropilWeight', 'var')
    NeuropilWeight = 0.65;
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numFramesBaseline'
                numFramesBaseline = varargin{index+1};
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

fprintf('Calculating trial-wise dF/F...');

%% Load ROI data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end
numROIs = numel(ROIdata.rois);

%% Determine number of frames to average for baseline
% if numFramesBaseline > ROIdata.DataInfo.numFramesBefore
%     numFramesBaseline = ROIdata.DataInfo.numFramesBefore;
% end

%% Calculate dF/F

% Initialize output
for rindex = 1:numROIs
    ROIdata.rois(rindex).dFoF = zeros(size(ROIdata.rois(rindex).data));
end

% Compute dF/F for each trial for each ROI
for rindex = 1:numROIs
    
    % Extract all trials for current stimulus
    data = ROIdata.rois(rindex).data;
    if NeuropilWeight % remove weighted Neuropil signal
        data = data - NeuropilWeight*ROIdata.rois(rindex).neuropil;
    end
    
    % Compute Fluorescence baseline for each trial
    baseline = median(data(:, ROIdata.DataInfo.numFramesBefore - numFramesBaseline + 1:ROIdata.DataInfo.numFramesBefore), 2);
    
    % Compute dF/F signal for each trial
    ROIdata.rois(rindex).dFoF = bsxfun(@rdivide, bsxfun(@minus, data, baseline), baseline);
end

fprintf('\tComplete\n');

%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end