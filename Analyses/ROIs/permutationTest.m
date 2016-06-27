function [dCoM,dSel,p,Actual] = permutationTest(ROIs, varargin)

numPerms = 10000;
distBetween = 1; % CoM
StimIDs = [];

ROIindex = [1,inf];
% TrialIndex = [1,inf];

saveOut = false;
saveFile = '';
UserData = [];

verbose = true;
directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
                index = index + 2;
%             case {'TrialIndex', 'Trials', 'trials'}
%                 TrialIndex = varargin{index+1};
%                 index = index + 2;
            case 'numPerms'
                numPerms = varargin{index+1};
                index = index + 2;
            case 'DistBtwn'
                distBetween = varargin{index+1};
                index = index + 2;
            case 'StimIDs'
                StimIDs = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
                index = index + 1;
            case {'save','Save'}
                saveOut = true;
                index = index + 1;
            case 'saveFile'
                saveFile = varargin{index+1};
                index = index + 2;
            case {'userdata','UserData'}
                UserData = varargin{index+1};
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

if ~exist('ROIs', 'var')
    [ROIs, p] = uigetfile({'*.rois'}, 'Select 2 corresponding ROI files', directory,'MultiSelect','on');
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p, ROIs);
end
if numel(ROIs) ~= 2
    error('Requires 2 ROI datasets: Full Pad data and its corresponding Single Whisker data');
end


%% Load data
if iscellstr(ROIs)
    ROIFiles = ROIs;
    for findex = 1:2
        load(ROIFiles{findex}, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata; clear ROIdata;
    end
end

% Determine file to save to
if saveOut && isempty(saveFile)
    if exist('ROIFiles','var')
        saveFile = strcat(strtok(ROIFiles{1},'.'),'.perm');
    else
        saveOut = false;
    end
end

% Compute mean stim
for findex = 1:2
    if ~isfield(ROIs{findex}.rois, 'stimMean')
        ROIs{findex} = computeTrialMean(ROIs{findex});
    end
end


%% Determine data to analyze

% Determine ROIs
if isrow(ROIindex)
    ROIindex = ROIindex';
end
if size(ROIindex,2) == 1
    ROIindex = repmat(ROIindex,1,2);
end
if ROIindex(end,1) == inf
    ROIindex = cat(1, ROIindex(1:end-1,:), repmat((ROIindex(1:end-1,1)+1:numel(ROIs{1}.rois))',1,2));
end
numROIs = size(ROIindex,1);

% % Determine trials
% if ~iscell(TrialIndex)
%     TrialIndex = {TrialIndex};
% end
% if numel(TrialIndex) == 1
%     TrialIndex = repmat(TrialIndex,2,1);
% end
% for findex = 1:2
%     if TrialIndex{findex}(end) == inf
%         TrialIndex{findex} = cat(2, TrialIndex{findex}(1:end-1), TrialIndex{findex}(1:end-1)+1:max(ROIs{findex}.DataInfo.TrialIndex));
%     end
%     TrialIndex{findex} = ismember(ROIs{findex}.DataInfo.TrialIndex', TrialIndex{findex});
% end

fprintf('Computing %d permutations for %d ROIs...\n',numPerms,numROIs);


%% Determine trials associated with each stimulus

% Determine stimuli
% if isempty(StimIDs)
    StimIDs = unique([ROIs{1}.DataInfo.StimID;ROIs{2}.DataInfo.StimID]);
% end
numStim = numel(StimIDs);

% % Determine trials for each stimulus
% Trials = cell(numStim,2);
% for findex = 1:2
%     for sindex = 1:numStim
%         Trials{sindex,findex} = ROIs{findex}.DataInfo.StimID==StimIDs(sindex); % determine all trials for current stimulus
%         Trials{sindex,findex}(~TrialIndex{findex}) = false;                    % keep only requested trials
%         Trials{sindex,findex} = find(Trials{sindex,findex});
%     end
% end
% numTrialsPerStim = cellfun(@numel,Trials);
% 
% % Offset and combine Trial indices
% Indices = cell(numStim,1);
% for sindex = 1:numStim
%     Indices{sindex} = [Trials{sindex,1};Trials{sindex,2}+numel(ROIs{1}.DataInfo.TrialIndex)];
% end
% numTrialsPerStim = cumsum(numTrialsPerStim,2);


%% Compute metrics off of permutations (trials same for all ROIs)

% % Create data matrix of trial means
% Data = [ROIs{1}.rois(ROIindex(:,1)).stimMean];
% Data = cat(1,Data,[ROIs{2}.rois(ROIindex(:,2)).stimMean]);
% if any(isnan(Data))
%     error('Data contains NaN(s)');
% end
% StimIDs = cat(1,ROIs{1}.DataInfo.StimID,ROIs{2}.DataInfo.StimID); %for debugging
% 
% % Initialize output
% dCoM = nan(numROIs,numPerms+1);
% dSel = nan(numROIs,numPerms+1);
% 
% % Cycle through permutations
% pfH = parfor_progress(numPerms+1);
% parfor pindex = 1:numPerms+1
%     
%     % Generate tuning curves (ignore control position)
%     Curves = nan(numROIs,numStim-1,2);
%     for sindex = 2:numStim
%         if pindex ~= 1
%             TI = Indices{sindex}(randperm(numTrialsPerStim(sindex,2)));                 % permute all trials for current stimulus
%         else
%             TI = Indices{sindex};                                                       % compute actual value
%         end
%         Curves(:,sindex-1,1) = mean(Data(TI(1:numTrialsPerStim(sindex,1)),:),1);        % average N trial means for each ROI
%         Curves(:,sindex-1,2) = mean(Data(TI(numTrialsPerStim(sindex,1)+1:end),:),1);    % average remaining trial means for second set of curves        
%     end
%     
%     % Compute metrics
%     dCoM(:,pindex) = diff([computeCenterOfMass(Curves(:,:,1),1,distBetween),computeCenterOfMass(Curves(:,:,2),1,distBetween)],[],2);
%     dSel(:,pindex) = diff([computeVectorSelectivity(Curves(:,:,1),[],1),computeVectorSelectivity(Curves(:,:,2),[],1)],[],2);
%     
%     parfor_progress(pfH); % update status
% end
% parfor_progress(pfH,0);
% 
% Actual = [dCoM(:,1),dSel(:,1)];
% dCoM(:,1) = [];
% dSel(:,1) = [];


%% Compute metrics off of permutations (w/ outliers removed)

% Gather Data
Data = [ROIs{1}.rois(ROIindex(:,1)).Raw,ROIs{2}.rois(ROIindex(:,1)).Raw];
Data = reshape(Data,numStim,numROIs,2);
numTrialsPerStim = cellfun(@numel,Data);
numTrialsPerStim = cumsum(numTrialsPerStim,3);
for rindex = 1:numROIs
    for sindex = 1:numStim
        Data{sindex,rindex,1} = [Data{sindex,rindex,1};Data{sindex,rindex,2}];
    end
end
Data(:,:,2) = [];

% Convert from cell to matrix (improves parfor)
num = max(max(numTrialsPerStim(:,:,2)));
temp = nan(numStim,numROIs,num);
for rindex = 1:numROIs
    for sindex = 1:numStim
       temp(sindex,rindex,:) = [Data{sindex,rindex};nan(num-numTrialsPerStim(sindex,rindex,2),1)];
    end
end
Data = temp;

% Initialize output
dCoM = nan(numROIs,numPerms+1);
dSel = nan(numROIs,numPerms+1);

% Cycle through permutations
pfH = parfor_progress(numPerms+1);
parfor pindex = 1:numPerms+1
    
    % Generate tuning curves (ignore control position)
    Curves = nan(numROIs,numStim-1,2);
    for sindex = 2:numStim
        for rindex = 1:numROIs
            if pindex ~= 1
                TI = randperm(numTrialsPerStim(sindex,rindex,2));
            else
                TI = 1:numTrialsPerStim(sindex,rindex,2);
            end
            Curves(rindex,sindex-1,1) = mean(Data(sindex,rindex,TI(1:numTrialsPerStim(sindex,rindex,1))));
            Curves(rindex,sindex-1,2) = mean(Data(sindex,rindex,TI(numTrialsPerStim(sindex,rindex,1)+1:numTrialsPerStim(sindex,rindex,2))));
        end
    end
    
    % Compute metrics
    dCoM(:,pindex) = diff([computeCenterOfMass(Curves(:,:,1),1,distBetween),computeCenterOfMass(Curves(:,:,2),1,distBetween)],[],2);
    dSel(:,pindex) = diff([computeVectorSelectivity(Curves(:,:,1),[],1),computeVectorSelectivity(Curves(:,:,2),[],1)],[],2);
    
    parfor_progress(pfH);
end
parfor_progress(pfH,0);

Actual = [dCoM(:,1),dSel(:,1)];
dCoM(:,1) = [];
dSel(:,1) = [];


%% Generate p-values
p = nan(numROIs,2);
p(:,1) = sum(bsxfun(@gt,abs(dCoM),abs(Actual(:,1))),2)/numPerms;
p(:,2) = sum(bsxfun(@gt,abs(dSel),abs(Actual(:,2))),2)/numPerms;


%% Save outputs
if saveOut
    if ~exist(saveFile,'file')
        save(saveFile,'dCoM','dSel','p','Actual','ROIindex','-v7.3');
    else
        save(saveFile,'dCoM','dSel','p','Actual','ROIindex','-append');
    end
    fprintf('Saved permutation results to: %s\n',saveFile);
end


