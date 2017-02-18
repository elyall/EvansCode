function [p,rho,Speed,theta,p_corr,CoM,p_tuned] = permutationRun(ROIdata, targetTrials, varargin)

numPerms = 5000;
StimIDs = [];       % vector of indices of stimuli to analyze
ROIindex = [];     
TrialIndex = [];    % indices of trials to analyze
RunSpeed = [];      % mean run speed for each trial
outliers = [];      % logical matrix specifying which trials are outliers for a given ROI (totalTrials x totalROIs)
pixelSize = [1,1];  % microns per pixel ([y,x]) -> corrects for aspect ratio
nSamples = [];

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
            case 'nSamples'
                nSamples = varargin{index+1};
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
            case {'saveFile','SaveFile'}
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
if ~iscell(targetTrials)
    targetTrials = {targetTrials};
end
numTargets = numel(targetTrials);
temp = false(totalTrials,numStim,numTargets);
for tindex = 1:numTargets
    temp(:,:,tindex) = bsxfun(@and, ismember(ROIdata.DataInfo.TrialIndex',targetTrials{tindex}), Stims);
end
targetTrials = temp;
if isempty(nSamples)
    nSamples = sum(targetTrials(:,:,1));
end


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
StimMeans = repmat(permute([ROIdata.rois(ROIindex).stimMean],[1,3,2]),1,numStim,1); % collect data (nTrials x nStim x nROIs)
StimMeans(~repmat(Stims,1,1,numROIs)) = NaN;                                        % remove data for other stimuli
StimMeans(repmat(permute(outliers,[1,3,2]),1,numStim,1)) = NaN;                     % remove outliers


%% Compute metrics off of permutations

% Initialize output
CoM = nan(numROIs,numPerms+numTargets);
TrialIndex = false(totalTrials,numPerms+numTargets);
p_tuned = nan(numROIs,numPerms+numTargets);
Max = nan(numROIs,numPerms+numTargets);

% Cycle through permutations computing map
blankTrials = false(totalTrials,numStim);
fprintf('Computing CoM for %d ROIs (%d permutations)...\n',numROIs,numPerms+numTargets);
pfH = parfor_progress(numPerms+numTargets);
parfor pindex = 1:numPerms+numTargets
    
    % Determine sample of trials for current permutation
    targetTrials; % parfor requires this variable be sent to all workers, otherwise returns indexing error
    if pindex <= numTargets
        currentTrials = targetTrials(:,:,pindex);
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
    
    % Determine peak of tuning curve
    Max(:,pindex) = max(Curves,[],2);
    
    % Compute center of mass 
    CoM(:,pindex) = computeCenterOfMass(Curves,positions,distBetween);
    
    % Compute p values for whether neurons are significantly tuned
    for rindex = 1:numROIs
        data = StimMeans(:,:,rindex);
        dict = arrayfun(@(x,y) y*ones(x,1),sum(currentTrials),1:numStim,'UniformOutput',false);
        p_tuned(rindex,pindex) = anovan(data(currentTrials), cat(1,dict{:}), 'model','full','display','off');
    end
    
    parfor_progress(pfH); % update status
end
parfor_progress(pfH,0);


%% Determine runspeed mean and variance for each permutation
if ~isempty(RunSpeed)
    SpeedMean = nan(1,numPerms+numTargets);
    SpeedStD = nan(1,numPerms+numTargets);
    parfor pindex = 1:numPerms+numTargets
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
index = ~isnan(CoM) & Max>.2 & p_tuned<.05;

% Initialize outputs
theta = nan(1,numPerms+numTargets);
rho = nan(1,numPerms+numTargets);
p_corr = nan(1,numPerms+numTargets);

fprintf('Computing map correlation for %d ROIs (%d permutations)...\n',numROIs,numPerms+numTargets);
pfH = parfor_progress(numPerms+numTargets);
parfor pindex = 1:numPerms+numTargets
    
    % Compute map correlation
    [Projection,theta(pindex)] = projectOntoCoMAxis(Centroids(index(:,pindex),:), CoM(index(:,pindex),pindex));
    [rho(pindex),p_corr(pindex)] = corr(Projection(:,1),CoM(index(:,pindex),pindex));
    
    parfor_progress(pfH); % update status
end
parfor_progress(pfH,0);

% Generate p-value
temp = mean(rho);
rho = rho - temp; % mean subtract distribution
p = nan(1,numTargets);
for tindex = 1:numTargets
    p(tindex) = sum(abs(rho(numTargets+1:end))>abs(rho(tindex)))/numPerms; % determine percent of permutations beyond target
    if p(tindex)>.0001
        fprintf('p%d=%.4f\t',tindex,p(tindex));
    elseif p(tindex)==0
        fprintf('p%d=%d',tindex,p(tindex));
    else
        fprintf('p%d=%e\t',tindex,p(tindex));
    end
end
fprintf('\n');
rho = rho + temp; % add mean back


%% Save outputs
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile,'file')
        save(saveFile,'p','rho','Speed','theta','p_corr','CoM','ROIindex','TrialIndex','numTargets','-v7.3');
    else
        save(saveFile,'p','rho','Speed','theta','p_corr','CoM','ROIindex','TrialIndex','numTargets','-append');
    end
    fprintf('Saved permutation results to: %s\n',saveFile);
end


%% Display plots
if verbose
    cmap = lines(numTargets);
    figure;
    subplot(1,3,1);
    plot(Speed.mean(numTargets+1:end),rho(numTargets+1:end),'k.'); hold on;
    for tindex = 1:numTargets
        plot(Speed.mean(tindex),rho(tindex),'*','Color',cmap(tindex,:));
    end
    legend([{'null dist.'};strcat('target',cellstr(num2str((1:numTargets)')))]);
    xlabel('Mean Run Speed (deg/s)');
    ylabel('Rho');
%     figure;
    subplot(1,3,2);
    plot(Speed.std(numTargets+1:end),rho(numTargets+1:end),'k.'); hold on;
    for tindex = 1:numTargets
        plot(Speed.std(tindex),rho(tindex),'*','Color',cmap(tindex,:));
    end
    legend([{'null dist.'};strcat('target',cellstr(num2str((1:numTargets)')))]);
    xlabel('STD Run Speed (deg/s)');
    ylabel('Rho');
%     figure;
    subplot(1,3,3);
    [N,edges] = histcounts(theta(numTargets+1:end),'Normalization','probability');
    stairs(edges([1,1:end,end]),[0,100*N([1:end,end]),0],'k');  hold on; % take edges down to 0 and fill out last bin
    YLim = get(gca,'YLim');
    for tindex = 1:numTargets
        plot([theta(tindex),theta(tindex)],YLim,'--','Color',cmap(tindex,:));
    end
    ylabel('Percent')
    xlabel('Angle of Axis (degrees)');
end
