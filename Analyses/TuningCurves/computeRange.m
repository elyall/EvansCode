function Range = computeRange(Curves, positions)


if ~exist('positions', 'var') || isempty(positions)
    positions = 2;
end


%% Compute range
numROIs = size(Curves,1);
if size(positions,1)==1
    positions = repmat(positions,numROIs,1);
end

Range = nan(numROIs,1);
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
        
    % Compute center of mass (preference)
    Range(rindex) = range(current);

end