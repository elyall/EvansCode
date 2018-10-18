function [ROIdata, Curves, outliers, p_tuned] = computeTuningCurve(ROIdata, ROIindex, TrialIndex, varargin)

FitTuningCurves = false; % gaussian fit
DetermineOutliers = true;
ControlID = 0; % StimID of control trials, or '[]' if no control trial
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
                FitTuningCurves = ~FitTuningCurves;
                index = index + 1;
            case {'Outliers','DetermineOutliers','outliers'}
                DetermineOutliers = ~DetermineOutliers;
                index = index + 1;
            case 'ControlID'
                ControlID = varargin{index+1};
                index = index + 2;
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

% Determine outliers for running trials of each stimulus
if DetermineOutliers
    fprintf('Determining outliers...');
    outliers = false(numel(TrialIndex),numel(ROIdata.rois));
    for rindex = ROIindex
        outliers(TrialIndex,rindex) = determineOutliers(ROIdata.rois(rindex).stimMean(TrialIndex),'GroupID',ROIdata.DataInfo.StimID(TrialIndex),'type','medianRule');
    end
    fprintf('\tComplete\n');
end

%% Determine stimuli info
if isempty(StimIDs)
    StimIDs = unique(ROIdata.DataInfo.StimID)';
elseif iscolumn(StimIDs)
    StimIDs = StimIDs';
end
numStimuli = numel(StimIDs);

% Determine order to cycle through stimuli (necessary for t-test to control trials)
if ~isempty(ControlID)
    controlindex = find(StimIDs==ControlID);                                    % locate control ID
    if ~isempty(controlindex)
        StimIDs = StimIDs([controlindex,setdiff(1:numStimuli,controlindex)]);   % reorder so control trials are analyzed first
    end
    firstStim = 2; % first index for significant tuning analysis (ANOVA)
else
    firstStim = 1; % first index for significant tuning analysis (ANOVA)
end

% Determine trials per stim
Trials = repmat(ROIdata.DataInfo.StimID,1,numStimuli)==repmat(StimIDs,numel(ROIdata.DataInfo.StimID),1); % determine in what trials each stimulus occurred
Trials = bsxfun(@and,Trials,TrialIndex); % keep only requested trials


%% Calculate average response for each stimulus
fprintf('Computing tuning curves...');

% Initialize output
[ROIdata.rois(ROIindex).curve] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).StdError] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).CI95] = deal(nan(2,numStimuli));
[ROIdata.rois(ROIindex).nTrials] = deal(nan(1, numStimuli));
[ROIdata.rois(ROIindex).Raw] = deal(cell(numStimuli, 1));
if ~isempty(ControlID)
    [ROIdata.rois(ROIindex).PValue] = deal(nan(1, numStimuli));
    [ROIdata.rois(ROIindex).PValueCorrected] = deal(nan(1, numStimuli-1));
end
[ROIdata.rois(ROIindex).TunedPValue] = deal(nan);

% Calculate tuning
for rindex = ROIindex
    ROIdata.rois(rindex).stimindex = StimIDs;
    for sindex = 1:numStimuli
                
        % Select data for current stimulus (ignore outliers)
        StimulusDFoF = ROIdata.rois(rindex).stimMean(Trials(:,sindex) & ~outliers(:,rindex));
        StimulusDFoF(isnan(StimulusDFoF)) = []; % remove nan trials (ROI didn't exist in trial either due to motion or bad merge across datasets)
        
        % Save tuning curves
        ROIdata.rois(rindex).curve(sindex) = mean(StimulusDFoF);                             % mean evoked dF/F over all trials for current stimulus
        ROIdata.rois(rindex).StdError(sindex) = std(StimulusDFoF)/sqrt(numel(StimulusDFoF)); % standard error of the mean
        ROIdata.rois(rindex).CI95(:,sindex) = bootci(10000,{@mean,StimulusDFoF},'type','bca'); % bootstrapped confidence intervals of the mean
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
    
    % Compute whether ROI is "tuned"
    if ~any(isnan(ROIdata.rois(rindex).curve))
        N = cellfun(@numel, ROIdata.rois(rindex).Raw(firstStim:end));
        dict = repelem(1:numel(N),N)';
        ROIdata.rois(rindex).TunedPValue = anovan(cat(1,ROIdata.rois(rindex).Raw{firstStim:end}), dict, 'model','full','display','off');
    else
        ROIdata.rois(rindex).TunedPValue = nan;
    end
        
end %ROIs
fprintf('\tComplete\n');


%% Fit tuning curves with gaussian
if FitTuningCurves
    fprintf('Fitting tuning curves with gaussian...');
    warning('off', 'curvefit:checkbounds:tooManyLowerBounds');
    for rindex = ROIindex
        if ~any(isnan(ROIdata.rois(rindex).curve))
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
        end
    end %ROIs
    fprintf('\tComplete\n');
end


%% Output extras
if nargout > 1
    Curves = reshape([ROIdata.rois(ROIindex).curve],numel(ROIdata.rois(1).curve),numel(ROIindex))';
end
if nargout > 3
    p_tuned = [ROIdata.rois(ROIindex).TunedPValue];
end


%% Save to file
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', 'outliers', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', 'outliers', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end
