function [TMI,x] = computeTMIbyPos(CurvesFull, CurvesSingle, positions, distBetween)

if ~exist('positions', 'var') || isempty(positions)
    positions = 2;
end

if ~exist('distBetween','var') || isempty(distBetween)
    distBetween = 1;
end


%% Compute preference
numROIs = size(CurvesFull,1);
if size(positions,1)==1
    positions = repmat(positions,numROIs,1);
end
if numel(distBetween)==1
    distBetween = repmat(distBetween,numROIs,1);
end

if ~iscell(CurvesFull)
    numStim = size(CurvesFull,2);
else
    numStim = max(cellfun(@numel,CurvesFull));
end

TMI = nan(numROIs,numStim);
for rindex = 1:numROIs
    
    % Select current ROI
    if ~iscell(CurvesFull)
        full = CurvesFull(rindex,:);
        single = CurvesSingle(rindex,:);
    else
        full = CurvesFull{rindex}';
        single = CurvesSingle{rindex}';
    end
    
    % Keep only region requested
    if ~iscell(positions) && size(positions,2)==1
        full = full(positions(rindex):end);
        full = full(positions(rindex):end);
    elseif iscell(positions)
        full = full(positions{rindex});
        full = full(positions{rindex});
    else
        full = full(positions(rindex,:));
        full = full(positions(rindex,:));
    end
    
    % Take magnitude of curve
    full = abs(full);
    full = abs(full);
    
    % Compute center of mass (preference)
    TMI(rindex,:) = x;

    
end