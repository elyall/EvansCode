function [Ratio, Center, Surround] = computeCenterVsSurround(Curves,PosCenter,PosSurround,Min)


%% Determine indexing
numROIs = size(Curves,1);

if size(PosCenter,1)==1
    PosCenter = repmat(PosCenter,numROIs,1);
end
if size(PosSurround,1)==1
    PosSurround = repmat(PosSurround,numROIs,1);
end
if numel(Min)==1
    Min = repmat(Min,numROIs,1);
end


%% Compute ratio

% Initialize output
Ratio = nan(numROIs,1);
Center = nan(numROIs,1);
Surround = nan(numROIs,1);

% Cycle through ROIs
for rindex = 1:numROIs
    
    % Select current ROI
    if ~iscell(Curves)
        current = Curves(rindex,:);
    else
        current = Curves{rindex}';
    end
    
    % shift bottom of curves to 0
    if isnumeric(Min) && ~isempty(Min);
        current = current - Min(rindex);
    elseif islogical(Min)
        current = current - min(current);
    end
    
    % Compute Center
    if iscell(PosCenter)
        Center(rindex) = mean(current(PosCenter{rindex}));
    else
        Center(rindex) = mean(current(PosCenter(rindex,1)));
    end
    
    % Compute Surround
    if iscell(PosSurround)
        Surround(rindex) = mean(current(PosSurround{rindex}));
    else
        Surround(rindex) = mean(current(PosSurround(rindex,1)));
    end
    
    % Compute Ratio
    Ratio(rindex) = Center(rindex)./Surround(rindex);

end