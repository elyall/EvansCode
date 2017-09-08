function s = lifetimeSparsity(Data,StimIndex)
% Data is ROI x trial or ROIdata or cell array of ROIdatas from multiple
% depths
% See: Froudarakis et al 2014 Nature; Vinje & Gallant 2008 Science

subtractMin = true;
verbose = false;

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

%% Compute lifetime sparsity of each ROI

% Compute mean across all trials for each stimulus
if exist('StimIndex','var') && ~isempty(StimIndex)
    [StimIDs,~,StimIndex] = unique(StimIndex);
    nS = numel(StimIDs);
    Data = arrayfun(@(x) mean(Data(:,StimIndex==x),2),1:nS,'UniformOutput',false);
    Data = cat(2,Data{:});
else
    nS = size(Data,2); % assume each column is it's own stimulus
end

% Min-subtract tuning curves
if subtractMin
    Data = bsxfun(@minus,Data,min(Data,[],2)); % subtract off minimum value
    % neg = any(Data<0,2);                                            % determine units that have negative responses
    % Data(neg,:) = bsxfun(@minus,Data(neg,:),min(Data(neg,:),[],2)); % subtract off minimum value for units that have responses below 0
end

% Compute lifetime sparseness for each ROI
s = (1-1/nS*sum(Data,2).^2./sum(Data.^2,2))/(1-1/nS);


%% Display results ordered from largest to smallest sparseness
if verbose
    [~,order] = sort(s,'descend');
    figure;
    imagesc(zscore(Data(order,:),[],2));
    xlabel('Stimulus'); ylabel('ROI');
    hCb = colorbar; ylabel(hCb,'z-score');
end
