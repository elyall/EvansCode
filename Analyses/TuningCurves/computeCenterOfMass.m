function CoM = computeCenterOfMass(Curves, positions, distBetween)

if ~exist('positions', 'var') || isempty(positions)
    positions = 2;
end

if ~exist('distBetween','var') || isempty(distBetween)
    distBetween = 1;
end


%% Compute preference
numROIs = size(Curves,1);
if size(positions,1)==1
    positions = repmat(positions,numROIs,1);
end
if numel(distBetween)==1
    distBetween = repmat(distBetween,numROIs,1);
end

CoM = nan(numROIs,1);
for rindex = 1:numROIs
    
    % Select current ROI
    if ~iscell(Curves)
        current = Curves(rindex,:);
    else
        current = Curves{rindex}';
    end
    
    % Keep only region requested
    if ~iscell(positions) && size(positions,2)==1
        current = current(positions(rindex):end);
    elseif iscell(positions)
        current = current(positions{rindex});
    else
        current = current(positions(rindex,:));
    end
    
    % Take magnitude of curve
    current = abs(current);
    
    % Compute center of mass (preference)
    index = distBetween(rindex):distBetween(rindex):numel(current)*distBetween(rindex); 
    CoM(rindex) = sum(index.*current)/sum(current);

    
end