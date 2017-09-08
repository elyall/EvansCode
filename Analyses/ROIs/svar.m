function SV = svar(Data,StimIndex)
% Data is ROI x trial or ROIdata or cell array of ROIdatas from multiple
% depths
% See: Froudarakis et al 2014 Nature

%% Parse input arguments
if isstruct(Data)
    Data = {Data};
end
if iscell(Data)
    first = Data{1}.DataInfo.numFramesBefore+1; % first stimulus frame
    last = Data{1}.DataInfo.numFramesBefore+mode(Data{1}.DataInfo.numStimFrames); % last stimulus frame
    StimIndex = Data{1}.DataInfo.StimID;
    temp = Data;
    Data = [];
    for findex = 1:numel(temp)
        Data = cat(3,Data,cat(3,temp{findex}.rois(:).dFoF));
    end
    Data = squeeze(mean(Data(:,first:last,:),2))'; % average over time and reshape to be ROI x trial
end

%% Compute spherical variance for the population

% Compute mean across all trials for each stimulus
if exist('StimIndex','var') && ~isempty(StimIndex)
    [StimIDs,~,StimIndex] = unique(StimIndex);
    nS = numel(StimIDs);
    Data = arrayfun(@(x) mean(Data(:,StimIndex==x),2),1:nS,'UniformOutput',false);
    Data = cat(2,Data{:});
else
    nS = size(Data,2); % assume each column is it's own stimulus
end

% Compute circular variance
s = nan(size(Data,1),nS);
for sindex = 1:nS
    s(:,sindex) = Data(:,sindex)/norm(Data(:,sindex)); % perform element-wise operation for each stimulus
end
SV = 1 - norm(1/nS*sum(s,2));

