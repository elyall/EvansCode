function outliers = determineOutliers(data, varargin)

numSTDsOutlier = 4; % inf if no threshold
minNumTrials = 5; % -inf if no minimum (will error when gets to 2 trials)
GroupID = [];

saveOut = false;
saveFile = '';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numSTDsOutlier'
                numSTDsOutlier = varargin{index+1};
                index = index + 2;
            case 'minNumTrials'
                minNumTrials = varargin{index+1};
                index = index + 2;
            case 'GroupID'
                GroupID = varargin{index+1};
                index = index + 2;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end


%% Initialize output (and break)
outliers = false(numel(data),1);
if numSTDsOutlier == inf
    return
end


%% Determine groups
if isempty(GroupID)
    GroupID = ones(numel(data),1);
end
IDs = unique(GroupID);


%% Determine outliers
for gindex = 1:nnz(~isnan(IDs))
    index = GroupID==IDs(gindex);
    outliers(index) = findOutliers(data(index), minNumTrials, numSTDsOutlier);
end


%% Save output
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile,'file')
        save(saveFile,'outliers','-mat','-v7.3');
    else
        save(saveFile,'outliers','-append');
    end
end


%% Main function
function outliers = findOutliers(data, minNumTrials, numSTDsOutlier)
N = numel(data);
outliers = false(N,1);

while true
    
    % Determine if the minimum number of trials is met
    n = nnz(~outliers);
    if n <= minNumTrials
        break
    end
    
    % Calculate z-score for each value relative to the rest of the
    % distribution
    Val = nan(N,1);
    for t = find(~outliers)'
        index = true(N,1);
        index(t) = false;
        index(outliers) = false;
        
        mu = nanmean(data(index));          % determine mean without current trial
        sigma = nanstd(data(index));        % determine std without current trial
        Val(t) = abs(data(t)-mu)/sigma;     % calculate zscore of current trial relative to other data
    end
    
    % Determine if any value is an outlier
    if any(Val > numSTDsOutlier)            % at least one outlier exists
        [~,furthestIndex] = max(Val);       % determine largest outlier
        outliers(furthestIndex) = true;     % remove largest outlier
    else
        break
    end
    
end

