function KFoldIndices = determineKFolds(TrialIndex)

numKFolds = 10;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numKFolds'
                numKFolds = varargin{index+1};
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


%% Determine parameters
if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:numel(ROIdata.DataInfo.StimID)];
end

    
%% Determine indices for k-means cross validation
KFoldIndices = nan(numel(TrialIndex), 2);
for sindex = 1:numStims
    
    % Determine in what trials current stimulus was shown
    currentIndices = find(StimIndex==sindex);
    numTrials = numel(currentIndices);
    
    % Determine number of trials that should exist in each fold
    numPerFold = diff(round(linspace(1,numTrials+1,numKFolds+1)));
    % numPerFold = numPerFold(randperm(numKFolds)); % randomize order of fold sizes (not necessary I believe)
    if any(numPerFold==0)
        warning('Requesting more folds than observations for a given stimulus -> some folds will not have a certain stimulus type.');
    end
    
    % Distribute the trials to the folds
    for kindex = 1:numKFolds
        temp = randsample(numTrials-sum(numPerFold(1:kindex-1)), numPerFold(kindex)); % select random trial(s) for current fold
        KFoldIndices(currentIndices(temp),:) = repmat([kindex, sindex], numPerFold(kindex), 1);
        currentIndices(temp) = [];
    end
    
end
    
