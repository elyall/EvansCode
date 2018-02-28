function s = lifetimeSparsity(Data,StimIndex,varargin)
% Data is ROI x trial or ROIdata or cell array of ROIdatas from multiple
% depths
% See: Froudarakis et al 2014 Nature; Vinje & Gallant 2008 Science

subtractMin = true;
fix = ''; % '', 'rectify', 'or 'truncate'
verbose = false;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'subtractMin'
                subtractMin = varargin{index+1};
                index = index + 2;
            case 'fix'
                fix = varargin{index+1};
                index = index + 2;
            case {'Verbose','verbose'}
                verbose = ~verbose;
                index = index + 1;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

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
end
nS = sum(~isnan(Data),2); % assume each column is it's own trial-averaged stimulus

% Min-subtract tuning curves
if isequal(subtractMin,true)
    Data = bsxfun(@minus,Data,min(Data,[],2)); % subtract off minimum value
elseif ~isempty(subtractMin)
    if numel(subtractMin)==1 && size(Data,1)~=1
        subtractMin = repmat(subtractMin,size(Data,1),1);
    end
    Data = bsxfun(@minus,Data,subtractMin); % subtract off input value(s)
end

% Deal with negative values
switch fix
    case 'rectify'
        Data = abs(Data); % make all negative values positive
    case 'truncate'
        Data(Data<0) = 0; % turn all negative values to 0
end

% Compute lifetime sparseness for each ROI
s = (1-1./nS.*nansum(Data,2).^2./nansum(Data.^2,2))./(1-1./nS);


%% Display results ordered from largest to smallest sparseness
if verbose
    [~,order] = sort(s,'descend');
    figure;
    imagesc(zscore(Data(order,:),[],2));
    xlabel('Stimulus'); ylabel('ROI');
    hCb = colorbar; ylabel(hCb,'z-score');
    
    %Data = bsxfun(@minus, Data, min(Data,[],2));
    %Data = bsxfun(@rdivide, Data, max(Data,[],2));
    %Data = bsxfun(@plus, Data, (size(Data,1)-.5:-1:.5)');
    figure;
    h = plot(Data(order,:)');
    set(h, {'Color'}, num2cell(flip(parula(numel(order))), 2));
end
