function Sel = computeVectorSelectivity(Curves, Min, firstpos)

if ~exist('min', 'var')
    Min = []; %empty means compute minimum individually
end

if ~exist('firstpos', 'var') || isempty(firstpos)
    firstpos = 2; % account for control position
end

%% Compute selectivity
numROIs = size(Curves,1);
Sel = nan(numROIs,1);

for rindex = 1:numROIs
    
    if ~iscell(Curves)
        current = Curves(rindex,:);
    else
        current = Curves{rindex}';
    end
    
    % Keep only selected region
    current = current(:,firstpos:end);
    
    % shift bottom of curves to 0
    if isempty(Min);
        current = current - min(current);
    elseif numel(Min)==1
        current = current - Min;
    else
        current = current - Min(rindex);
    end
    
    % compute vector selectivity
    Sel(rindex) = 1 - (sqrt(sum(abs(current).^2))/max(current)-1)/(sqrt(numel(current))-1); % (sqrt(|x1|^2+|x2|^2)/max(x) - 1 / sqrt(numel(x))) - 1
    
end
