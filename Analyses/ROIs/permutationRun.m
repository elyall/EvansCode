function [p,rho,Speed,theta,p_corr,CoM] = permutationRun(ROIdata, targetTrials, varargin)

numPerms = 500;
StimIDs = [];       % vector of indices of stimuli to analyze
ROIindex = [];     
TrialIndex = [];    % indices of trials to analyze
RunSpeed = [];      % mean run speed for each trial
outliers = [];      % logical matrix specifying which trials are outliers for a given ROI (totalTrials x totalROIs)
pixelSize = [1,1];  % microns per pixel ([y,x]) -> corrects for aspect ratio

% CoM calculation
positions = [];
distBetween = 1;

verbose = true;

saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numPerms'
                numPerms = varargin{index+1};
                index = index + 2;
            case 'StimIDs'
                StimIDs = varargin{index+1};
                index = index + 2;
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'TrialIndex', 'Trials', 'trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'outliers'
                outliers = varargin{index+1};
                index = index + 2;
            case 'RunSpeed'
                RunSpeed = varargin{index+1};
                index = index + 2;
            case 'pixelSize'
                pixelSize = varargin{index+1};
                index = index + 2;
            case 'positions'
                positions = varargin{index+1};
                index = index + 2;
            case {'DistBtwn','distBetween'}
                distBetween = varargin{index+1};
                index = index + 2;
            case {'save','Save'}
                saveOut = true;
                index = index + 1;
            case 'saveFile'
                saveFile = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
                index = index + 1;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('ROIdata', 'var')
    [ROIdata, p] = uigetfile({'*.rois'}, 'Select 2 corresponding ROI files', directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p, ROIdata);
end


%% Load data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', 'outliers', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = strcat(strtok(ROIFile,'.'),'.permrun');
    end
end

% Compute mean stim
if ~isfield(ROIdata.rois, 'stimMean')
    ROIdata = computeTrialMean(ROIdata);
end


%% Determine data to analyze

% Determine trials to use
totalTrials = numel(ROIdata.DataInfo.TrialIndex);
if isempty(TrialIndex) % use all trials
    TrialIndex = true(totalTrials,1);
else
    TrialIndex = ismember(ROIdata.DataInfo.TrialIndex',TrialIndex); % convert to logical vector
end

% Determine stimuli to use
if isempty(StimIDs)
    StimIDs = unique(ROIdata.DataInfo.StimID);
    StimIDs(1) = []; % remove control stimulus
end
numStim = numel(StimIDs);

% Determine trial index of each stimulus\
Stims = bsxfun(@and, TrialIndex, bsxfun(@eq, ROIdata.DataInfo.StimID, StimIDs'));  % logical vector indexing trials for each stimulus
StimIndex = cell(1,numStim);
for sindex = 1:numStim
    StimIndex{sindex} = find(Stims(:,sindex));
end

% Determine number of trials to pull for each stimulus
targetTrials = ismember(ROIdata.DataInfo.TrialIndex',targetTrials); % logical vector indexing target trials
targetTrials = bsxfun(@and, targetTrials, Stims);                   % logical vector numTrials x numStim indexing target trials for each stimulus
nSamples = sum(targetTrials);


%% Determine ROIs to analyze and gather their responses
if isempty(ROIindex)
    ROIindex = 1:numel(ROIdata.rois);
end
numROIs = numel(ROIindex);

% Determine outliers
if exist('ROIFile','var') && ~isempty(outliers) % outliers were loaded
    outliers = outliers(:,ROIindex);            % remove unwanted ROIs
elseif isempty(outliers)                        % no outliers thrown out
    outliers = false(totalTrials,numROIs);
elseif isequal(outliers, true)                  % compute outliers
    outliers = false(totalTrials,numROIs);
    for rindex = 1:numROIs
        outliers(TrialIndex,rindex) = determineOutliers(ROIdata.rois(ROIindex(rindex)).stimMean(TrialIndex),'GroupID',ROIdata.DataInfo.StimID(TrialIndex),'type','medianRule');
    end
end

% Gather mean activity for each trial of every ROI
StimMeans = repmat(permute([ROIdata.rois(ROIindex).stimMean],[1,3,2]),1,numStim,1); % collect data
StimMeans(~repmat(Stims,1,1,numROIs)) = NaN;                                        % remove data for other stimuli
StimMeans(repmat(permute(outliers,[1,3,2]),1,numStim,1)) = NaN;                     % remove outliers


%% Compute metrics off of permutations

% Initialize output
CoM = nan(numROIs,numPerms+1);
TrialIndex = false(totalTrials,numPerms+1);

% Cycle through permutations computing map
blankTrials = false(totalTrials,numStim);
fprintf('Computing CoM for %d ROIs (%d permutations)...\n',numROIs,numPerms);
pfH = parfor_progress(numPerms+1);
for pindex = 1:numPerms+1
    
    % Determine trials to sample
    if pindex == 1
        currentTrials = targetTrials;
    else
        currentTrials = blankTrials;
        for sindex = 1:numStim
            currentTrials(randsample(StimIndex{sindex},nSamples(sindex)),sindex) = true;
        end
    end
    TrialIndex(:,pindex) = any(currentTrials,2);
    
    % Generate tuning curves (ignore control position)
    currentData = StimMeans;
    currentData(~repmat(currentTrials,1,1,numROIs)) = NaN;
    Curves = squeeze(nanmean(currentData))';
    Curves(all(Curves<=.2,2),:) = NaN; % remove not driven neurons
    
    % Compute center of mass 
    CoM(:,pindex) = computeCenterOfMass(Curves,positions,distBetween);
         
    parfor_progress(pfH); % update status
end
parfor_progress(pfH,0);


%% Determine runspeed mean and variance for each permutation
if ~isempty(RunSpeed)
    SpeedMean = nan(1,numPerms+1);
    SpeedStD = nan(1,numPerms+1);
    parfor pindex = 1:numPerms+1
        SpeedMean(pindex) = mean(RunSpeed(TrialIndex(:,pindex)));
        SpeedStD(pindex) = std(RunSpeed(TrialIndex(:,pindex)));
    end
    Speed.mean = SpeedMean;
    Speed.std = SpeedStD;
    clear SpeedMean SpeedStD
else
    Speed = [];
end


%% Compute correlations of map

% Determine centroids
Centroids = cat(1,ROIdata.rois(ROIindex).centroid);     % gather centroids (units in pixels)
Centroids = bsxfun(@times,Centroids,pixelSize([2,1]));  % convert units to um (corrects for aspect ratio)
index = ~isnan(CoM);

% Initialize outputs
theta = nan(1,numPerms+1);
rho = nan(1,numPerms+1);
p_corr = nan(1,numPerms+1);

fprintf('Computing map correlation for %d ROIs (%d permutations)...\n',numROIs,numPerms);
pfH = parfor_progress(numPerms+1);
for pindex = 1:numPerms+1
    
    % Compute map correlation
    [Projection,theta(pindex)] = projectOntoCoMAxis(Centroids(index(:,pindex),:), CoM(index(:,pindex),pindex));
    [rho(pindex),p_corr(pindex)] = corr(Projection(:,1),CoM(index(:,pindex),pindex));
    
    parfor_progress(pfH); % update status
end
parfor_progress(pfH,0);

% Generate p-value
temp = mean(rho);
rho = rho - temp;
p = sum(abs(rho(2:end))>abs(rho(1)))/numPerms;
rho = rho + temp;
if p>.0001
    fprintf('p = %.4f\n',p);
else
    fprintf('p = %e\n',p);
end

%% Save outputs
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile,'file')
        save(saveFile,'p','rho','Speed','theta','p_corr','CoM','ROIindex','TrialIndex','-v7.3');
    else
        save(saveFile,'p','rho','Speed','theta','p_corr','CoM','ROIindex','TrialIndex','-append');
    end
    fprintf('Saved permutation results to: %s\n',saveFile);
end


%% Display plots
if verbose
    figure;
    subplot(1,3,1);
    plot(Speed.mean(2:end),rho(2:end),'k.'); hold on;
    plot(Speed.mean(1),rho(1),'g*');
    xlabel('Mean Run Speed (deg/s)');
    ylabel('Rho');
%     figure;
    subplot(1,3,2);
    plot(Speed.std(2:end),rho(2:end),'k.'); hold on;
    plot(Speed.std(1),rho(1),'g*');
    xlabel('STD Run Speed (deg/s)');
    ylabel('Rho');
%     figure;
    subplot(1,3,3);
    histogram(theta,'Normalization','probability');
    ylabel('Percent')
    xlabel('Angle of Axis (degrees)');
end
