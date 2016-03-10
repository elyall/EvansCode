function CoM = computeCoM2(Curves)

xDist = 1;
yDist = 1;

%% Compute preference
[H,W,numROIs] = size(Curves);
CoM = nan(numROIs,2);

for rindex = 1:numROIs
    
    % Select current ROI
    if ~iscell(Curves)
        current = Curves(:,:,rindex);
    else
        current = Curves{rindex}';
    end
    
    % Take magnitude of curve
    current = current>0;
    
    % Compute center of mass (preference)
    [X,Y] = meshgrid(xDist:xDist:W*xDist, yDist:yDist:H*yDist);
    CoM(rindex,:) = [sum(X(:).*current(:))/sum(current(:)), sum(Y(:).*current(:))/sum(current(:))];

    
end