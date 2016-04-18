function [ROIdata, goodTrials, Curves] = computeTuningCurve(ROIdata, ROIindex, TrialIndex, varargin)


FitTuningCurves = false; % gaussian fit
ControlID = 0; % StimID of control trials, or '[]' if no control trial
outlierweight = 3; % # of std dev to ignore
StimIDs = [];

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
            case 'StimIDs'
                StimIDs = varargin{index+1};
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

% Compute trial means
if ~isfield(ROIdata.rois, 'stimMean')
    ROIdata = computeTrialMean(ROIdata);
end


%% Determine data to analyze

if TrialIndex(end) == inf
    TrialIndex = cat(2, TrialIndex(1:end-1), TrialIndex(1:end-1)+1:max(ROIdata.DataInfo.TrialIndex));
end
TrialIndex = ismember(ROIdata.DataInfo.TrialIndex', TrialIndex);
% if isrow(TrialIndex)
%     TrialIndex = TrialIndex';
% end

if ischar(ROIindex) && strcmp(ROIindex, 'all')
    ROIindex = [1, inf];
elseif iscolumn(ROIindex)
    ROIindex = ROIindex';
end
if ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:numel(ROIdata.rois));
end


%% Determine stimuli info
if isempty(StimIDs)
    StimIDs = unique(ROIdata.DataInfo.StimID(TrialIndex));
end
numStimuli = numel(StimIDs);

% Determine order to cycle through stimuli (necessary for t-test to control trials)
if ~isempty(ControlID)
    controlindex = find(StimIDs==ControlID);                                    % locate control ID
    if ~isempty(controlindex)
        StimIDs = StimIDs([controlindex,setdiff(1:numStimuli,controlindex)]);   % reorder so control trials are analyzed first
    end
end

% Determine trials per stim
Trials = repmat(ROIdata.DataInfo.StimID,1,numStimuli)==repmat(StimIDs,numel(ROIdata.DataInfo.StimID),1); % determine in what trials each stimulus occurred
Trials = bsxfun(@and,Trials,TrialIndex); % keep only requested trials
TrialIndex = cell(numStimuli,1);
% for sindex=1:numStimuli
%     TrialIndex{sindex} = find(Trials(:,sindex));
% end


%% Calculate average response for each stimulus
fprintf('Computing tuning curves...');

% Initialize output
[ROIdata.rois(ROIindex).curve] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).StdError] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).nTrials] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).Raw] = deal(cell(numStimuli, 1));
if ~isempty(ControlID)
    [ROIdata.rois(ROIindex).PValue] = deal(nan(1, numStimuli));
    [ROIdata.rois(ROIindex).PValueCorrected] = deal(nan(1, numStimuli-1));
end

% Calculate tuning
% goodTrials = repmat(TrialIndex,1,numROIs);
for rindex = ROIindex
    ROIdata.rois(rindex).stimindex = StimIDs;
    for sindex = 1:numStimuli
                
        % Select data for current stimulus
        StimulusDFoF = ROIdata.rois(rindex).stimMean(Trials(:,sindex));
        
        % Remove outliers
        while true
            numTrials = numel(StimulusDFoF);
            Val = nan(numTrials,1);
            for tindex = 1:numTrials
                mu = mean(StimulusDFoF(setdiff(1:numTrials,tindex)));               % determine mean without current trial
                sigma = std(StimulusDFoF(setdiff(1:numTrials,tindex)));             % determine std without current trial
                Val(tindex) = abs(StimulusDFoF(tindex)-mu)/sigma;                   % calculate zscore of current trial relative to other data
            end
            if any(Val > outlierweight)                                             % at least one outlier exists
                [~,furthestIndex] = max(Val);                                       % determine largest outlier
                StimulusDFoF(furthestIndex) = [];                                   % remove largest outlier
                % goodTrials{sindex,rindex}(furthestIndex) = [];
            else
                break
            end
        end
        % ROIdata.rois(rindex).stimMean(setdiff(TrialIndex{sindex},goodTrials{sindex,rindex})) = nan;
                
        % Save tuning curves
        ROIdata.rois(rindex).curve(sindex) = nanmean(StimulusDFoF);                 % evoked dF/F over all trials for current stimulus
        ROIdata.rois(rindex).StdError(sindex) = std(StimulusDFoF)/sqrt(numel(StimulusDFoF)); % standard error for stimulus
        ROIdata.rois(rindex).Raw{sindex} = StimulusDFoF;
        ROIdata.rois(rindex).nTrials(sindex) = numel(StimulusDFoF);
        
        % Perform t-test to control trials
        if ~isempty(ControlID)
            if sindex==1 %control trials
                ControlDFoF = StimulusDFoF;
            else
                [~, ROIdata.rois(rindex).PValue(sindex)] = ttest2(...
                    ControlDFoF,....
                    StimulusDFoF);
            end
        end
        
    end %stimuli
    
    % Correct for multiple comparisons
    if ~isempty(ControlID)
        [~,~,ROIdata.rois(rindex).PValueCorrected] = fdr_bh(ROIdata.rois(rindex).PValue);
    end
        
end %ROIs
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
        ROIdata.rois(rindex).FWHM = ROIdata.rois(rindex).Coeff(3)*2*sqrt(2*log(2)); % FWHM = 2*sqrt(2*ln(2))*c (c is std dev of curve)
        ROIdata.rois(rindex).rsquare = ROIdata.rois(rindex).GoFit.rsquare;
        
    end %ROIs
    fprintf('\tComplete\n');
end


%% Output curves
if nargout > 2
    Curves = reshape([ROIdata.rois(ROIindex).curve],numel(ROIdata.rois(1).curve),numel(ROIindex))';
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
