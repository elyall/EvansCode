function Corr = computeCurveCorrelation(Curves,positions,DataIndex)

if ~exist('positions','var')
    positions = [];
end
if ~iscell(positions)
    positions = {positions};
end

if ~exist('DataIndex','var') || isempty(DataIndex)
    DataIndex = ones(size(Curves,1),1);
end

IDs = unique(DataIndex);
numMice = numel(IDs);

if numel(positions)==1 && numMice>1
    positions = repmat(positions,numMice,1);
end

%% Pearson's
Corr = [];
for Mindex = 1:numMice
    
    % Pull out tuning curves
    if iscell(Curves)
        current = cell2mat(Curves(DataIndex==IDs(Mindex))');
    else
        current = Curves(DataIndex==IDs(Mindex),:)';
    end
    
    % Keep only requested positions
    if ~isempty(positions{Mindex})
        current = current(positions{Mindex},:);
    end
    
    % compute pearson's r
    temp = corrcoef(current);
    Corr = [Corr;temp(tril(true(size(temp)),-1))];
    
end


%% Kendall's