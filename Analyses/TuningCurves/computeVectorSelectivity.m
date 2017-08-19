function Sel = computeVectorSelectivity(Curves, Min, positions)

if ~exist('min', 'var')
    Min = []; %empty means compute minimum individually
end

if ~exist('positions', 'var') || isempty(positions)
    positions = 2;
end

%% Compute selectivity
numROIs = size(Curves,1);
if size(positions,1)==1
    positions = repmat(positions,numROIs,1);
end
if numel(Min)==1
    Min = repmat(Min,numROIs,1);
end

Sel = nan(numROIs,1);
for rindex = 1:numROIs
    
    if ~iscell(Curves)
        current = Curves(rindex,:);
    else
        current = Curves{rindex}';
    end
    
    % Keep only region requested
    if isnumeric(positions) && size(positions,2)==1
        current = current(positions(rindex):end);
    elseif iscell(positions)
        current = current(positions{rindex});
    else
        current = current(positions(rindex,:));
    end
    
    % shift bottom of curves to 0
    if isempty(Min);
        current = current - min(current);
    else
        current = current - Min(rindex);
    end
    
    % compute vector selectivity
    Sel(rindex) = 1 - (sqrt(sum(abs(current).^2))/max(current)-1)/(sqrt(numel(current))-1); % 1 - (sqrt(|x1|^2+|x2|^2)/max(x) - 1) / (sqrt(numel(x)) - 1)
    
end
