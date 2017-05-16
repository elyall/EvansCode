function [ROIdata,Means] = computeTrialMean(ROIdata, varargin)

ROIindex = [1 inf];
FrameIndex = []; % default takes mean over whole stim period; set to 'VertPole' to analyze last 500ms of stim

saveOut = false;
saveFile = '';

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

% Compute dF/F
if ~isfield(ROIdata.rois, 'dFoF')
    ROIdata = computeDFoF(ROIdata);
end
numTrials = size(ROIdata.rois(1).dFoF,1);

% Fix legacy issue
if ~isfield(ROIdata,'Config')
    warning('ROIdata doesn''t have Config struct -> assuming FrameRate=15.49 & Depth=1');
    ROIdata.Config.Depth = 1;
    ROIdata.Config.FrameRate = 15.49;
end

%% Determine data to analyze
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end

if isempty(FrameIndex)
    FrameIndex = cell(numTrials,1);
    for tindex = 1:numTrials
        FrameIndex{tindex} = ROIdata.DataInfo.numFramesBefore+1:ROIdata.DataInfo.numFramesBefore+ROIdata.DataInfo.numStimFrames(tindex); % analyze whole stim period
    end
elseif ischar(FrameIndex) && strcmp(FrameIndex,'VertPole')
    FrameIndex = cell(numTrials,1);
    for tindex = 1:numTrials
        F = ROIdata.DataInfo.numFramesBefore+ROIdata.DataInfo.numStimFrames(tindex);
        FrameIndex{tindex} = F-floor(round(ROIdata.Config.FrameRate/2)/ROIdata.Config.Depth)+1:F; % analyze last 500ms
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
        ROIdata.rois(rindex).stimMean(tindex) = nanmean(ROIdata.rois(rindex).dFoF(tindex,FrameIndex{tindex}),2); % nan frames may exist due to motion artifacts
    end
end
fprintf('\tComplete\n');


%% Output means
if nargout > 1
    Means = reshape([ROIdata.rois(ROIindex).stimMean],numTrials,numel(ROIindex))';
end


%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end
