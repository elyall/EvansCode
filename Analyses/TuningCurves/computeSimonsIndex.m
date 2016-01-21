function SI = computeSimonsIndex(Curves, firstpos)

if ~exist('firstpos', 'var') || isempty(firstpos)
    firstpos = 1; % account for control position
end

if ~iscell(Curves)
    
    % shift bottom of curves to 0
    Curves = bsxfun(@minus, Curves(:,firstpos:end), min(Curves(:,firstpos:end), [], 2));
    
    % compute index
    SI = max(Curves, [], 2)./mean(Curves, 2);
    
else
   
    % shift bottom of curves to 0
    Curves = cellfun(@(x) x(firstpos:end)-min(x(firstpos:end)), Curves, 'UniformOutput', false);
    
    % compute index
    SI = cellfun(@(x) max(x)/mean(x), Curves);
    
end