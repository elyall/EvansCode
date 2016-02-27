function CoM = computeCenterOfMass(Curves, firstpos)

distBetween = 1;

if ~exist('firstpos', 'var') || isempty(firstpos)
    firstpos = 2;
end

%% Compute preference
numROIs = size(Curves,1);
CoM = nan(numROIs,1);

for rindex = 1:numROIs
    
    % Select current ROI
    if ~iscell(Curves)
        current = Curves(rindex,:);
    else
        current = Curves{rindex}';
    end
    
    % Keep only region requested
    if exist('firstpos','var')
        current = current(:,firstpos:end);
    end
    
    % Take magnitude of curve
    current = abs(current);
    
    % Compute center of mass (preference)
    index = distBetween:distBetween:numel(current)*distBetween; 
    CoM(rindex) = sum(index.*current)/sum(current);

    
end