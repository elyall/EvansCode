function outliers = determineOutliers(data, varargin)

type = 'medianRule'; % 'tukey','evan','medianRule',

% Evan's method
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
            case {'type','Type'}
                type = varargin{index+1};
                index = index + 2;
            case 'numSTDsOutlier'
                numSTDsOutlier = varargin{index+1};
                index = index + 2;
            case 'minNumTrials'
                minNumTrials = varargin{index+1};
                index = index + 2;
            case 'GroupID'
                GroupID = varargin{index+1};
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


%% Initialize output (and break)
outliers = false(numel(data),1);
if strcmp(type,'evan') && numSTDsOutlier == inf
    return
end


%% Determine groups
if isempty(GroupID)
    GroupID = ones(numel(data),1);
end
IDs = unique(GroupID);


%% Determine outliers
for gindex = 1:nnz(~isnan(IDs))
    index = GroupID==IDs(gindex); % if nan exists, UNIQUE() places it at end
    switch type
        case 'evan'
            outliers(index) = findOutliers(data(index), minNumTrials, numSTDsOutlier);
        case 'tukey'
            outliers(index) = tukeyMethod(data(index));
        case 'medianRule'
            outliers(index) = medianRule(data(index));
    end
end


%% Save output
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile,'file')
        save(saveFile,'outliers','-mat','-v7.3');
    else
        save(saveFile,'outliers','-append');
    end
end


%% Main functions
% Tukey method
function outliers = tukeyMethod(data)
% q = quantile(data,[.25,.75]); % same thing
q = prctile(data,[25,75]);
IQR = diff(q);
edges = [q(1)-1.5*IQR, q(2)+1.5*IQR];
outliers = data < edges(1) | data > edges(2);

% Median rule
function outliers = medianRule(data)
m = median(data);
IQR = diff(prctile(data,[25,75]));
edges = [m-2.3*IQR, m+2.3*IQR];
outliers = data < edges(1) | data > edges(2);

% Evan method
function outliers = findOutliers(data, minNumTrials, numSTDsOutlier)
if iscolumn(data)
    data = data';
end
N = numel(data);

squaremat = repmat(data',1,N); % create square, diagonal matrix
squaremat(1:N+1:end) = nan;    % remove diagonal to exclude own value from relative zscore

outliers = false(N,1);
while true
    
    % Determine if the minimum number of trials is met
    if nnz(~outliers) <= minNumTrials
        break
    end
    
    % Calculate z-score for each value relative to the rest of the data
    ZScore = abs(data-nanmean(squaremat))./nanstd(squaremat);
    
    % Determine if any value is an outlier
    if any(ZScore > numSTDsOutlier)
        furthestIndex = find(Zscore==max(ZScore)); % determine largest outlier(s)
        outliers(furthestIndex) = true;            % record largest outlier
        squaremat(furthestIndex,:) = NaN;          % remove outlier from inclusion in other distributions
        squaremat(:,furthestIndex) = NaN;          % remove outlier as an option
    else
        break
    end
    
end

