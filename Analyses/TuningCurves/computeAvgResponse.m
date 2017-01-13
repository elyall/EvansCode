function [AvgResponse, StimIDs, outliers] = computeAvgResponse(ROIdata, ROIindex, TrialIndex, varargin)

StimIDs = [];

saveOut = false;
saveFile = '';

directory = cd;

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'StimIDs'
                StimIDs = varargin{index+1};
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
    [ROIdata, p] = uigetfile({'*.rois;*.mat'}, 'Choose ROI file', directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if ~exist('ROIindex','var') || isempty(ROIindex)
    ROIindex = [1 inf];
end

if ~exist('TrialIndex', 'var') || isempty(TrialIndex)
    TrialIndex = [1 inf];
elseif islogical(TrialIndex)
    TrialIndex = find(TrialIndex);
end


%% Load data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
end

% Compute trial means
if ~isfield(ROIdata.rois, 'stimMean')
    ROIdata = computeTrialMean(ROIdata);
end
numFrames = size(ROIdata.rois(1).dFoF,2);


%% Determine data to analyze
if TrialIndex(end) == inf
    TrialIndex = cat(2, TrialIndex(1:end-1), TrialIndex(end-1)+1:max(ROIdata.DataInfo.TrialIndex));
end
TrialIndex = ismember(ROIdata.DataInfo.TrialIndex', TrialIndex);
% if isrow(TrialIndex)
%     TrialIndex = TrialIndex';
% end

% Determine ROIs
if ischar(ROIindex) && strcmp(ROIindex, 'all')
    ROIindex = [1, inf];
elseif iscolumn(ROIindex)
    ROIindex = ROIindex';
end
if ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:numel(ROIdata.rois));
end
numROIs = numel(ROIindex);

% Determine outliers for running trials of each stimulus
fprintf('Determining outliers...');
outliers = false(numel(TrialIndex),numROIs);
for rindex = 1:numROIs
    outliers(TrialIndex,rindex) = determineOutliers(ROIdata.rois(ROIindex(rindex)).stimMean(TrialIndex),'GroupID',ROIdata.DataInfo.StimID(TrialIndex),'type','medianRule');
end
fprintf('\tComplete\n');


%% Determine stimuli info
allStimIDs = unique(ROIdata.DataInfo.StimID)';
if isempty(StimIDs)
    StimIDs = allStimIDs;
end
if ~iscell(StimIDs)
    StimIDs = {StimIDs};
end
if numel(StimIDs)==1
    StimIDs = repmat(StimIDs,numROIs,1);
end
numStimuli = max(cellfun(@numel,StimIDs));

% Determine trials per stim
Trials = repmat(ROIdata.DataInfo.StimID,1,numel(allStimIDs))==repmat(allStimIDs,numel(ROIdata.DataInfo.StimID),1); % determine in what trials each stimulus occurred
Trials = bsxfun(@and,Trials,TrialIndex); % keep only requested trials


%% Calculate average response for each stimulus
fprintf('Computing average response...');

% Initialize output
AvgResponse = nan(numROIs,numFrames,numStimuli);

% Calculate tuning
for rindex = 1:numROIs
    for sindex = 1:numel(StimIDs{rindex})
        
        % Select data for current stimulus (ignore outliers)
        StimulusDFoF = ROIdata.rois(ROIindex(rindex)).dFoF(Trials(:,StimIDs{rindex}(sindex)==allStimIDs) & ~outliers(:,rindex),:);
        StimulusDFoF(any(isnan(StimulusDFoF),2),:) = []; % remove nan trials (ROI didn't exist in trial either due to motion or bad merge across datasets)
        
        % Average trials together
        AvgResponse(rindex,:,sindex) = mean(StimulusDFoF);
        
    end %stimuli       
end %ROIs
fprintf('\tComplete\n');


%% Save to file
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile, 'file')
        save(saveFile, 'AvgResponse', '-mat', '-v7.3');
    else
        save(saveFile, 'AvgResponse', '-mat', '-append');
    end
    fprintf('\AvgResponse saved to: %s\n', saveFile);
end
