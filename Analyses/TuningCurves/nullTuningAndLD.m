function [Curves, SE, p, LD, CI, order] = nullTuningAndLD(Data, StimID, StimLog, varargin)

TrialIndex = [1,inf];
numIters = 100;
numBoot = 500;

saveOut = false;
saveFile = '';

%% Parse input data
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'TrialIndex'
                TrialIndex = varargin{index+1};
                index = index + 2;
            case {'numIters','N','num','iters'}
                numIters = varargin{index+1};
                index = index + 2;
            case 'numBoot'
                numBoot = varargin{index+1};
                index = index + 2;
            case {'Save','save','saveOut'}
                saveOut = varargin{index+1};
                index = index + 2;
            case 'saveFile'
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
        
%% Load in data
if ischar(Data)
    saveFile = [Data(1:end-4),'null'];
    load(Data,'ROIdata','-mat');
    Data = ROIdata;
end
if isstruct(Data)
    Data = gatherROIdata(Data,'stimMean');
end

if ischar(StimID)
    ExpFile = StimID;
    load(ExpFile,'TrialInfo','-mat');
    StimID = TrialInfo.StimID;
end

if ~exist('StimLog','var') && exist('ExpFile','var')
    load(ExpFile,'Experiment','-mat');
    StimLog = Experiment.stim.stim;
elseif ischar(StimLog)
    load(StimLog,'Experiment','-mat');
    StimLog = Experiment.stim.stim;
end

%% Keep only selected trials
if isinf(TrialIndex(end))
    TrialIndex = [TrialIndex(1:end-1),TrialIndex(end-1)+1:size(Data,2)];
end
Data = Data(:,TrialIndex);
StimID = StimID(TrialIndex);

%% Compute tuning curves and linear differences
numStims = size(StimLog,1);
numMW = nnz(sum(StimLog,2)>1);
[numROIs,numTrials] = size(Data);
Curves = nan(numROIs,numStims,numIters);
SE = nan(numROIs,numStims,numIters); %SE = nan(numROIs,2,numStims,numIters);
p = nan(numROIs,numStims,numIters);
LD = nan(numROIs,numMW,numIters);
CI = nan(numROIs,2,numMW,numIters);
order = nan(numTrials,numIters);
for ind = 1:numIters
    
    % Randomize stim labels
    order(:,ind) = randperm(numTrials);
    current = StimID(order(:,ind));
    
    % Compute tuning
    [Curves(:,:,ind), SE(:,:,ind), temp, Raw] = computeTuningCurve2(Data, current, 'verbose', false);
    p(:,:,ind) = temp.driven_p_corr;
    
    % Compute linear difference
    [LD(:,:,ind), CI(:,:,:,ind)] = LinDiff(Raw,StimLog,'N',numBoot);
    
end


%% Save ouput
if saveOut && ~isempty(saveFile)
    if exist(saveFile,'file')
        save(saveFile,'C','SE','p','pLD','pCI','ord','-append');
    else
        save(saveFile,'C','SE','p','pLD','pCI','ord','-v7.3');
    end
end

