function TMI = computeTMI(Curves, region)

numROIs = size(Curves,1);
TMI = nan(numROIs,1);

for rindex = 1:numROIs
    
    % Select current ROI
    if ~iscell(Curves)
        current = Curves(rindex,:,:);
    else
        current = cat(3, Curves{rindex,1}',Curves{rindex,2}');
    end
    
    % Keep only region requested
    if exist('region','var')
        if ~iscell(region)
            if size(region,1) == 1
                current = current(:,region,:);
            else
                current = current(:,region(rindex,:),:);
            end
        else
            if numel(region) == 1
                current = current(:,region{1},:);
            else
                current = current(:,region{rindex},:);
            end
        end
    end
    
    % Shift bottom of curves to 0
    current = current - min(current(:));
    
    % Take mean value across all positions input
    val = squeeze(mean(current,2));
    
    % Compute Trimming Modulation Index
    TMI(rindex) = diff(val)/sum(val);
    
end