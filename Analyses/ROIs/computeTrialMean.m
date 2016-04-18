function ROIdata = computeTrialMean(ROIdata, varargin)

saveOut = false;
saveFile = '';

ROIindex = [1 inf];
FrameIndex = [];

directory = cd;

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FrameIndex'
                FrameIndex = varargin{index+1};
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

if ~exist('ROIdata','var') || isempty(ROIdata)
    [ROIdata, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

fprintf('Calculating mean activity per trial...');


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

% Compute dF/F
if ~isfield(ROIdata.rois, 'dFoF')
    ROIdata = computeDFoF(ROIdata, 0.65);
end
numTrials = size(ROIdata.rois(1).dFoF,1);


%% Determine data to analyze
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end

if isempty(FrameIndex)
    FrameIndex = cell(numTrials,1);
    for tindex = 1:numTrials
        FrameIndex{tindex} = ROIdata.DataInfo.numFramesBefore+1:ROIdata.DataInfo.numFramesBefore+ROIdata.DataInfo.numStimFrames(tindex);
    end
elseif ~iscell(FrameIndex)
    FrameIndex = {FrameIndex};
end
if numel(FrameIndex) == 1
    FrameIndex = repmat(FrameIndex,numTrials,1);
end


%% Calculate mean per trial

% Initialize output
[ROIdata.rois(ROIindex).stimMean] = deal(nan(numTrials,1));

% Compute mean for each trial
for rindex = ROIindex
    for tindex = 1:numTrials
        ROIdata.rois(rindex).stimMean(tindex) = nanmean(ROIdata.rois(rindex).dFoF(tindex,FrameIndex{tindex}),2);
    end
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