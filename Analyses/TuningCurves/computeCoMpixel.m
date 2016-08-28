function CoM = computeCoMpixel(FoVStimAvg, positions, distBetween)

if ~exist('positions', 'var') || isempty(positions)
    positions = 2;
end

if ~exist('distBetween','var') || isempty(distBetween)
    distBetween = 1;
end

if numel(positions)==1
    FoVStimAvg = FoVStimAvg(:,:,positions:end);
elseif numel(positions)>1
    FoVStimAvg = FoVStimAvg(:,:,positions);
end

%% Compute preference
[H,W,S] = size(FoVStimAvg);
CoM = nan(H,W);
    
% Take magnitude of curve
FoVStimAvg = abs(FoVStimAvg);

% Compute center of mass (preference)
index = repmat(distBetween:distBetween:S*distBetween,H,1,W);
index = permute(index,[1,3,2]);
CoM = sum(index.*FoVStimAvg,3)./sum(FoVStimAvg,3);

