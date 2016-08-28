function [CoM, Sel, numToChoose] = permuteData(ROIdata, varargin)

type = 'multiple'; % 'split' or 'multiple'
perc = .5; % 'multiple' only
numPerms = 10000; % 'multiple' only

ROIindex = [1 inf];
TrialIndex = [1 inf];

PWCZ = 2;           % CoM
distBetween = 1;    % CoM

directory = cd;

% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'TrialIndex', 'Trials', 'trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case 'perc'
                perc = varargin{index+1};
                index = index + 2;
            case 'DistBtwn'
                distBetween = varargin{index+1};
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

if ~exist('ROIdata', 'var')
    [ROIdata, p] = uigetfile({'*.rois;*.mat'}, 'Select ROI file', directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p, ROIdata);
end


%% Load data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
end

% Compute mean stim
if ~isfield(ROIdata.rois, 'stimMean')
    ROIdata = computeTrialMean(ROIdata);
end


%% Determine data to analyze
if ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:numel(ROIdata.rois));
end
numROIs = numel(ROIindex);

if TrialIndex(end) == inf
    TrialIndex = cat(2, TrialIndex(1:end-1), TrialIndex(1:end-1)+1:numel(ROIdata.DataInfo.StimID));
end
numTrials = numel(TrialIndex);


%% Determine permutations

% Determine trials for each stimulus
StimIDs = unique(ROIdata.DataInfo.StimID(ismember(ROIdata.DataInfo.TrialIndex,TrialIndex)));
numStim = numel(StimIDs);
Trials = cell(numStim,1);
numToChoose = nan(numStim,1);
for sindex = 1:numStim
    Trials{sindex} = find(ROIdata.DataInfo.StimID==StimIDs(sindex));
    Trials{sindex}(~ismember(ROIdata.DataInfo.TrialIndex(Trials{sindex}),TrialIndex)) = [];
    numToChoose(sindex) = round(perc*numel(Trials{sindex}));
end

% Determine trials for each permutation
Permutations = cell(1,numStim);
switch type
    case 'split' 
        
        % Only use 2 permutations: first half and second half
        numPerms = 2;
        for sindex = 1:numStim
            num = floor(numel(Trials{sindex})/2);
            Permutations{sindex}(1,:) = Trials{sindex}(1:num);          %first half
            Permutations{sindex}(2,:) = Trials{sindex}(end-num+1:end);  %second half
        end
        numToChoose = [];
        
    case 'multiple'
        
        % Create permutations
        for sindex = 1:numStim
            Permutations{sindex} = nan(numPerms,numToChoose(sindex));
            for pindex = 1:numPerms
                Permutations{sindex}(pindex,:) = datasample(Trials{sindex}, numToChoose(sindex), 1, 'Replace', false);
            end
        end
        
end %type


%% Compute new data
Data = [ROIdata.rois(ROIindex).stimMean];

CoM = nan(numROIs,numPerms);
Sel = nan(numROIs,numPerms);

parfor_progress(numPerms);
parfor pindex = 1:numPerms
    Curves = nan(numROIs,numStim);
    for sindex = 1:numStim
        Curves(:,sindex) = mean(Data(Permutations{sindex}(pindex,:),:),1);
    end
    CoM(:,pindex) = computeCenterOfMass(Curves, PWCZ, distBetween);
    Sel(:,pindex) = computeVectorSelectivity(Curves, [], 2);
    parfor_progress;
end
parfor_progress(0);


