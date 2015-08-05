function [ROIdata, badTrials] = computeTuningCurve(ROIdata, ROIindex, TrialIndex, varargin)


FitTuningCurves = false; % gaussian fit
ControlID = 0; % StimID of control trials, or '[]' if no control trial
outlierweight = 3; % # of std dev to ignore
StimFrames = [];

saveOut = false;
saveFile = '';

directory = cd;

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Fit'
                FitTuningCurves = true;
                index = index + 1;
            case 'ControlID'
                ControlID = varargin{index+1};
                index = index + 2;
            case 'StimFrames'
                StimFrames = varargin{index+1};
                index = index + 2;
            case 'outlierweight'
                outlierweight = varargin{index+1};
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
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end

% Compute dF/F
if ~isfield(ROIdata.rois, 'dFoF')
    ROIdata = computeDFoF(ROIdata, 0.65);
end


%% Determine trials to analyze
if TrialIndex(end) == inf
    TrialIndex = cat(2, TrialIndex(1:end-1), TrialIndex(1:end-1)+1:numel(ROIdata.DataInfo.StimID));
end


%% Determine ROIs to analyze
if ischar(ROIindex) && strcmp(ROIindex, 'all')
    ROIindex = [1, inf];
end
if ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:numel(ROIdata.rois));
end

%% Determine stimuli info
StimIDs = unique(ROIdata.DataInfo.StimID(TrialIndex));
numStimuli = numel(StimIDs);

% Determine order to cycle through stimuli (necessary for t-test to control trials)
if ~isempty(ControlID)
    controlindex = find(StimIDs==ControlID); %locate control 'stimulus' in structure
    if ~isempty(controlindex)
        indexorder = cat(1, controlindex, find(StimIDs~=ControlID)); %run through the stimuli analyzing the control trials first
    else
        ControlID = [];
        indexorder = (1:numStimuli)';
    end
end

% Determine stimulus frames for each stimulus
if isempty(StimFrames)
    StimFrames = zeros(numStimuli, 2);
    StimFrames(:, 1) = ROIdata.DataInfo.numFramesBefore + 1;
    for sindex = 1:numStimuli
        StimFrames(sindex, 2) = StimFrames(1, 1) + mode(ROIdata.DataInfo.numStimFrames(ROIdata.DataInfo.StimID==StimIDs(sindex))) - 1;
    end
elseif isvector(StimFrames)
    StimFrames = repmat(StimFrames, numStimuli, 1);
end


%% Calculate average response for each stimulus
fprintf('Computing tuning curves...');

% Initialize output
[ROIdata.rois(ROIindex).curve] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).StdError] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).nTrials] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).Raw] = deal(cell(numStimuli, 1));
if ~isempty(ControlID)
    [ROIdata.rois(ROIindex).PValue] = deal(nan(1, numStimuli));
end

% Calculate tuning
badTrials = cell(numel(ROIdata.DataInfo.StimID), 1);
for rindex = ROIindex
    for sindex = indexorder'
        
        % Select data for current stimulus
        currentTrials = find(ROIdata.DataInfo.StimID==StimIDs(sindex));     % determine all trials for current stimulus
        currentTrials = currentTrials(ismember(currentTrials,TrialIndex));  % remove non-specified trials
        CaTraces = ROIdata.rois(rindex).dFoF(currentTrials,:);              % pull out data for these trials
        
        % Remove bad trials
        badCurrent = any(isnan(CaTraces), 2);
        if any(badCurrent)
            for bindex = currentTrials(badCurrent)
                badTrials{bindex} = [badTrials{bindex}, rindex];
            end
            CaTraces(badCurrent, :) = [];
        end
        
        % Compute response for each trial during the stimulus period
        StimulusDFoF = mean(CaTraces(:, StimFrames(sindex,1):StimFrames(sindex,2)), 2);
        
        % Remove outliers
        n = inf;
        while n ~= length(StimulusDFoF)
            n = length(StimulusDFoF);
            StimulusDFoF(abs(StimulusDFoF-mean(StimulusDFoF))>outlierweight*std(StimulusDFoF)) = [];
        end
        
        % Save tuning curves
        ROIdata.rois(rindex).curve(sindex) = mean(StimulusDFoF); % evoked dF/F over all trials for current stimulus
        ROIdata.rois(rindex).StdError(sindex) = std(StimulusDFoF)/sqrt(length(StimulusDFoF)); % standard error for stimulus
        ROIdata.rois(rindex).Raw{sindex} = StimulusDFoF;
        ROIdata.rois(rindex).nTrials(sindex) = numel(StimulusDFoF);
        
        % Perform t-test to control trials
        if ~isempty(ControlID)
            if StimIDs(sindex) == ControlID %control trials
                ControlDFoF = StimulusDFoF;
            end
            [~, ROIdata.rois(rindex).PValue(sindex)] = ttest2(...
                ControlDFoF,....
                StimulusDFoF);
        end
        
    end %stimuli
end %ROIs
if any(~cellfun(@isempty, badTrials))
    warning('Some trials couldn''t be included as they contained at least one NaN');
end
fprintf('\tComplete\n');


%% Fit tuning curves with gaussian
if FitTuningCurves
    fprintf('Fitting tuning curves with gaussian...');
    warning('off', 'curvefit:checkbounds:tooManyLowerBounds');
    for rindex = ROIindex
        
        % Ignore control trials for tuning curves
        if ~isempty(ControlID)
            minStimIndex = 2;
        else
            minStimIndex = 1;
        end
        
        % Determine DC component
        ROIdata.rois(rindex).min = min(ROIdata.rois(rindex).curve(minStimIndex:end));
        ROIdata.rois(rindex).max = max(ROIdata.rois(rindex).curve(minStimIndex:end));
        
        % Subtract off DC component in either direction
        upcurve = ROIdata.rois(rindex).curve(minStimIndex:end)-ROIdata.rois(rindex).min;
        downcurve = ROIdata.rois(rindex).curve(minStimIndex:end)-ROIdata.rois(rindex).max;

        % Fit upward and downward gaussian
        [UpFit, UpGoFit] = FitFull(upcurve);
        [DownFit, DownGoFit] = FitFull(downcurve);
        
        % Save fit paraemters
        if UpGoFit.rsquare >= DownGoFit.rsquare
            ROIdata.rois(rindex).FitDirection = 'Up';
            ROIdata.rois(rindex).Fit = UpFit;
            ROIdata.rois(rindex).GoFit = UpGoFit;
            ROIdata.rois(rindex).offset = ROIdata.rois(rindex).min;
        else
            ROIdata.rois(rindex).FitDirection = 'Down';
            ROIdata.rois(rindex).Fit = DownFit;
            ROIdata.rois(rindex).GoFit = DownGoFit;
            ROIdata.rois(rindex).offset = ROIdata.rois(rindex).max;
        end
        ROIdata.rois(rindex).Coeff = coeffvalues(ROIdata.rois(rindex).Fit);
        ROIdata.rois(rindex).ConfIntervals = confint(ROIdata.rois(rindex).Fit);
        ROIdata.rois(rindex).rsquare = ROIdata.rois(rindex).GoFit.rsquare;
        
    end %ROIs
    fprintf('\tComplete\n');
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
