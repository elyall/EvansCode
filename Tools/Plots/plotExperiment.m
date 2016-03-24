function Data = plotExperiment(Data, Neuropil, ExpStimFrames, varargin)

Type = 'line'; % 'line' or 'image'
ROIindex = [1 inf];
FrameIndex = [1 inf];

% Neuropil
baselinePrctile = 30;

% Smoothing
smoothType = ''; % '', 'moving', 'lowess', 'loess', 'sgolay', 'rlowess', or 'rloess'

hA = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Type'
                Type = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FrameIndex'
                FrameIndex = varargin{index+1};
                index = index + 2;
            case 'baselinePrctile'
                baselinePrctile = varargin{index+1};
                index = index + 2;
            case 'smoothType'
                smoothType = varargin{index+1};
                index = index + 2;
            case 'hA'
                hA = varargin{index+1};
                index = index + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Data','var') || isempty(Data)
    [Data, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file',directory);
    if isnumeric(Data)
        return
    end
    Data = fullfile(p,Data);
end


%% Load ROIs

% Load in ROIdata
if ischar(Data)
    ROIFile = Data;
    Data = load(ROIFile, 'ROIdata', '-mat');
    Data = Data.ROIdata;
end

% Convert ROIdata to matrices
if isstruct(Data)
    Neuropil = reshape([Data.rois(:).rawneuropil], numel(Data.rois(1).rawdata), numel(Data.rois));
    Data = reshape([Data.rois(:).rawdata], numel(Data.rois(1).rawdata), numel(Data.rois));
end
[totalFrames, totalROIs] = size(Data);

% Determine ROIs to compute weights for
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(1:end-1)+1:totalROIs];
end
numROIs = numel(ROIindex);
Data = Data(:,ROIindex);

% Determine frames to pull from
if FrameIndex(end) == inf
    FrameIndex = [FrameIndex(1:end-1), FrameIndex(1:end-1)+1:totalFrames];
end

% Format experiment data
if exist('ExpStimFrames','var') && ~isempty(ExpStimFrames)
    if isvector(ExpStimFrames)
        Stim = ExpStimFrames;
    else
        Stim = zeros(1,numel(FrameIndex));
        for tindex = 1:size(ExpStimFrames,1)
            Stim(ExpStimFrames(tindex,1):ExpStimFrames(tindex,2)) = 1;
        end
    end
    Stim = Stim(FrameIndex);
end


%% Subtract off neuropil
if exist('Neuropil', 'var') && ~isempty(Neuropil)
    NeuropilWeight = determineNeuropilWeight(Data, Neuropil(:,ROIindex));
    Data = Data - bsxfun(@times,NeuropilWeight',Neuropil(:,ROIindex));
end


%% Compute dFoF
if ~isempty(baselinePrctile)
    baseline = prctile(Data, baselinePrctile);
    Data = bsxfun(@rdivide, bsxfun(@minus, Data, baseline), baseline);
end


%% Smooth data
if ~isempty(smoothType)
    Data = smooth(Data, smoothType);
end


%% Plot activity
if isempty(hA)
    figure;
    hA = axes();
end

Data = Data(FrameIndex,:);
switch Type
    case 'image'
        imagesc(Data');
        
    case 'line'
        x = Data;
        if numROIs > 1 % scale
            x = bsxfun(@minus, x, min(x));
            x = bsxfun(@rdivide, x, max(x));
            x = bsxfun(@plus, x, 0.5:numel(ROIindex)-.5); % space out lines
        end
                
        % Plot stimuli
        if exist('Stim','var')
            area(any(Stim,1)*max(x(:)),'FaceColor',[.9,.9,.9],'EdgeColor',[.9,.9,.9]); hold on;
        end
        
        % Plot data
        plot(x);
                
end
