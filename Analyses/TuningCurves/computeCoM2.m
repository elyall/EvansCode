function CoM = computeCoM2(Curves)

xDist = 1;
yDist = 1;

%% Compute preference
[H,W,numROIs] = size(Curves);
CoM = nan(numROIs,2);

parfor rindex = 1:numROIs
    
    % Select current ROI
    if ~iscell(Curves)
        current = Curves(:,:,rindex);
    else
        current = Curves{rindex}';
    end
    
    % Take magnitude of curve
    current = abs(current);
    
    % Compute projections
    currentx = mean(current);
    currenty = mean(current,2)';
    
    % Compute center of mass (preference)
    CoM(rindex,:) = [sum(currentx.*(1:W))/sum(currentx),sum(currenty.*(1:H))/sum(currenty)];
%     [X,Y] = meshgrid(xDist:xDist:W*xDist, yDist:yDist:H*yDist);
%     CoM(rindex,:) = [sum(X(:).*current(:))/sum(current(:)), sum(Y(:).*current(:))/sum(current(:))];

    
end