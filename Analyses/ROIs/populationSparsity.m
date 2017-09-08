function s = populationSparsity(Data,StimIndex)
% Data is ROI x trial or ROIdata or cell array of ROIdatas from multiple
% depths
% See: Froudarakis et al 2014 Nature; Vinje & Gallant 2008 Science

subtractMin = true;

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

%% Compute population sparsity for each stimulus

% Compute mean across all trials for each stimulus
if exist('StimIndex','var') && ~isempty(StimIndex)
    [StimIDs,~,StimIndex] = unique(StimIndex);
    Data = arrayfun(@(x) mean(Data(:,StimIndex==x),2),1:numel(StimIDs),'UniformOutput',false);
    Data = cat(2,Data{:});
end

% Min-subtract tuning curves
if subtractMin
    Data = bsxfun(@minus,Data,min(Data,[],2)); % subtract off minimum value
    % neg = any(Data<0,2);                                            % determine units that have negative responses
    % Data(neg,:) = bsxfun(@minus,Data(neg,:),min(Data(neg,:),[],2)); % subtract off minimum value for units that have responses below 0
end

% Compute population sparseness for each stimulus
nR = size(Data,1);
s = (1-1/nR*sum(Data,1).^2./sum(Data.^2,1))/(1-1/nR);

