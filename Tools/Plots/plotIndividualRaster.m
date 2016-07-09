function [hA, CLim, CMap] = plotIndividualRaster(rois, ROIindex, DataInfo, varargin)

% Data to display
dataType = 'dFoF';      % 'dFoF' or 'raw'
sortType = 'stim';      % 'stim' or 'trial' (will not show averages)
showData = true;
showControl = false;    % nan means no control given
showAvgs = false;
TrialIndex = [1 inf];
stimOrder = [];

% Computation parameters
numSTDsOutlier = inf;   % 'inf' means no trials are thrown out

% Display parameters
zscoreby = 'none';      %'all' or 'trial' or 'none'
avgPixelHeight = 6;     % number of pixels of average plots height
frameRate = 15.45;      % frame rate of the data
xlab='Frames';          % default x-axis label
framelimit = 'all';     % xlim (indexes data points)
lineWidth = .5;
labelAxes = false;
XLabels = {};

% Color properties
showColorBar = false;
CLim = [];
CMap = {};

% Spike parameters
showSpikes = false;
spikeThresh = .18;

% Constant parameters
hA = []; %axes handle

directory = cd;

%% Parse input arguments
if ~exist('rois','var') || isempty(rois)
    [rois, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file',directory);
    if isnumeric(rois)
        return
    end
    rois = fullfile(p,rois);
end

if ~exist('ROIindex','var') || isempty(ROIindex)
    ROIindex = [1 inf];
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Trials','trials','TrialIndex'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'stimOrder'
                stimOrder = varargin{index+1};
                index = index + 2;
            case 'dataType'
                dataType = varargin{index+1};
                index = index + 2;
            case 'sortType'
                sortType = varargin{index+1};
                index = index + 2;
            case 'showData'
                showData = varargin{index+1};
                index = index + 2;
            case 'showAvgs'
                showAvgs = varargin{index+1};
                index = index + 2;
            case 'numSTDsOutlier'
                numSTDsOutlier = varargin{index+1};
                index = index + 2;
            case 'zscoreby'
                zscoreby = varargin{index+1};
                index = index + 2;
            case 'avgPixelHeight'
                avgPixelHeight = varargin{index+1};
                index = index + 2;
            case 'framelimit'
                framelimit = varargin{index+1};
                index = index + 2;
            case 'FrameRate'
                frameRate = varargin{index+1};
                xlab = 'Time (s)';
                index = index + 2;
            case 'lineWidth'
                lineWidth = varargin{index+1};
                index = index + 2;
            case 'axes'
                hA = varargin{index+1};
                index = index + 2;
            case 'showColorBar'
                showColorBar = true;
                index = index + 1;
            case 'CLim'
                CLim = varargin{index+1};
                index = index + 2;
            case 'CMap'
                CMap = varargin{index+1};
                index = index + 2;
            case 'showControl'
                showControl = true;
                index = index + 1;
            case 'labelAxes'
                labelAxes = varargin{index+1};
                index = index + 2;
            case 'XLabels'
                XLabels = varargin{index+1};
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

%% Load in data
if ischar(rois)
    ROIFile = rois;
    if saveToPDF && strcmp(saveFile, 'test.pdf')
        saveFile = strcat(ROIFile(1:end-3), 'pdf');
    end
    load(ROIFile, 'ROIdata', 'RunningSpeed', '-mat');
    rois = ROIdata.rois;
    DataInfo = ROIdata.DataInfo;
end
if isrow(DataInfo.TrialIndex)
    DataInfo.TrialIndex = DataInfo.TrialIndex';
end


%% Determine trials to display
if TrialIndex(end)==inf
    TrialIndex = cat(2, TrialIndex(1:end-1), TrialIndex(end-1)+1:size(rois(1).data, 1));
elseif islogical(TrialIndex)
    TrialIndex = find(TrialIndex);
end
TrialIndex = ismember(DataInfo.TrialIndex,TrialIndex);

StimIDs = unique(DataInfo.StimID(TrialIndex)); % determine stimuli present

% Determine whether to show control trials
if islogical(showControl) && ~showControl
    TrialIndex(DataInfo.StimID==StimIDs(1)) = false;
    StimIDs(1) = [];
end
numStims = numel(StimIDs);
numTrials = nnz(TrialIndex);

% Determine stimulus frames
StimIndex = DataInfo.StimID(TrialIndex);
StimFrames = nan(numTrials,2);
StimFrames(:,1) = DataInfo.numFramesBefore+1; % first frame of stimulus
% StimFrames(:,2) = DataInfo.numFramesBefore + DataInfo.numStimFrames(TrialIndex);
for sindex = 1:numStims
    StimFrames(StimIndex==StimIDs(sindex),2) = mode(DataInfo.numFramesBefore + DataInfo.numStimFrames(TrialIndex & DataInfo.StimID==StimIDs(sindex)));
end
if islogical(showControl) && showControl
    StimFrames(DataInfo.StimID(TrialIndex)==StimIDs(1),:) = nan; % do not show stim bars for control trials
end


%% Determine ROIs to plot
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(1:end-1)+1:numel(rois)];
end
numROIs = numel(ROIindex);


%% Determine color properties
if isempty(CLim)
    CLim = zeros(numROIs, 2);
elseif isvector(CLim)
    CLim = repmat(CLim, numROIs, 1);
end
if isempty(CMap)
    CMap = cell(numROIs, 1);
elseif isnumeric(CMap)
    CMap = repmat({CMap}, numROIs, 1);
elseif numel(CMap)==1
    CMap = repmat(CMap, numROIs, 1);
end


%% Create axes
if isempty(hA)
    for rindex = 1:numROIs
        figure();
        hA(rindex) = axes();
    end
end


%% Plot rasters
for rindex = 1:numROIs
        
    % Extract data
    switch dataType
        case 'dFoF'
            data = rois(ROIindex(rindex)).dFoF(TrialIndex, :);
        case 'raw'
            data = rois(ROIindex(rindex)).data(TrialIndex, :);
    end
    if showSpikes
        spikes = rois(ROIindex(rindex)).spikes(TrialIndex,:);
    end
    meanData = rois(ROIindex(rindex)).stimMean(TrialIndex);
    
    % Remove outliers for each stimulus
    if strcmp(sortType, 'stim') && numSTDsOutlier
        good = true(numTrials,1);
        for sindex = 1:numStims
            current = find(StimIndex==StimIDs(sindex));
            outliers = determineOutliers(meanData(current));
            good(current(outliers)) = false;
        end
        data(~good,:) = [];
        StimFrames(~good,:) = [];
        StimIndex(~good) = [];
        if showSpikes
            spikes(~good,:) = [];
        end
        numTrials = numTrials - nnz(~good);
    end
    
    % Determine number of trials per stimulus
    if strcmp(sortType, 'stim')
        blockLength = nan(numStims,1);
        for sindex = 1:numStims
            blockLength(sindex) = sum(StimIndex==StimIDs(sindex));
        end
    end
    
    % Determine sorting order
    switch sortType
        case 'stim'
            if isempty(stimOrder)
                [~, displayIndex] = sort(StimIndex);
            else
                displayIndex = nan(numTrials,1);
                stimInd = nan(numStims,1);
                ind = 0;
                for sindex=1:numStims
                    stimInd(sindex) = find(StimIDs==stimOrder(sindex));
                    displayIndex(ind+1:ind+blockLength(stimInd(sindex))) = find(StimIndex==stimOrder(sindex));
                    ind = ind + blockLength(StimIDs==stimOrder(sindex));
                end
                blockLength = blockLength(stimInd);
            end
            if isempty(XLabels)
                XLabels = cellstr(num2str((1:numStims)'));
                if islogical(showControl) && showControl
                    XLabels{1} = 'nc';
                end
            end

        case 'trial'
            displayIndex = 1:numTrials;
            showAvgs = false;
    end
    
    % Sort data
    data = data(displayIndex,:);
    StimFrames = StimFrames(displayIndex,:);
    StimIndex = StimIndex(displayIndex);
    if showSpikes
        spikes = spikes(displayIndex,:);
    end
        
%     if showAvgs
%         % Compute average stimuli
%         AvgStim = zeros(numStims*avgPixelHeight, size(data, 2));
%         AvgStimuliFrames = [];
%         for sindex = 1:numStims
%             AvgStim((sindex-1)*avgPixelHeight+1:(sindex-1)*avgPixelHeight+avgPixelHeight, :) = repmat(mean(data(TrialID==sindex, :), 1), avgPixelHeight, 1);
%             AvgStimuliFrames = cat(1, AvgStimuliFrames, repmat(mode(StimFrames(StimIndex==sindex, :), 1), avgPixelHeight, 1));
%         end
%         
%         % Add average stimuli to matrix to display
%         if showData
%             data = cat(1, data, AvgStim);
%             blockLength = cat(1, blockLength, repmat(avgPixelHeight, numStims,1));
%             StimFrames = cat(1, StimFrames, AvgStimuliFrames);
%             if ~isnan(ControlID)
%                 labels = cat(1, labels, [{'nc avg'}; cellstr([num2str(StimIDs(~ismember(StimIDs, ControlID))), repmat(' avg', numStims-1,1)])]);
%             else
%                 labels = cat(1, labels, cellstr([num2str(StimIDs), repmat(' avg', numStims,1)]));
%             end
%         else
%             data = AvgStim;
%             blockLength = repmat(avgPixelHeight, numStims,1);
%             StimFrames = AvgStimuliFrames;
%             if ~isnan(ControlID)
%                 labels = [{'nc'}, cellstr([num2str(StimIDs(~ismember(StimIDs, ControlID))), repmat(' avg', numStims-1,1)])];
%             else
%                 labels = cellstr([num2str(StimIDs), repmat(' avg', numStims,1)]);
%             end
%         end
%     end
    
    % Z-score data
    switch zscoreby
        case 'trial'
            data = zscore(data,0,2); %zscore each row individually
            colorbarTitle = 'z-score (per trial)';
        case 'all'
            [H,W] = size(data);
            data = reshape(zscore(data(:)),H,W);
            colorbarTitle = 'z-score (all)';
        case 'none'
            % do nothing
            switch dataType
                case 'raw'
                    colorbarTitle = 'Fluorescence (AU)';
                case 'dFoF'
                    colorbarTitle = 'dF/F';
            end
    end
    
    % Crop data
    if isnumeric(framelimit) && numel(framelimit)==2
        xshift = -framelimit(1)+1;
        data = data(:,framelimit(1):framelimit(2));
        if showSpikes
            spikes = spikes(:,framelimit(1):framelimit(2));
        end
    else
        xshift = 0;
    end
    
    % Determine color attritbutes
    if range(CLim(rindex,:)) == 0
        CLim(rindex,:) = [min(data(:)), max(data(:))];
    end
    if isempty(CMap{rindex})
        if ~strcmp(zscoreby, 'none') || strcmp(dataType, 'dFoF')
            CMap{rindex} = b2r(CLim(rindex,1),CLim(rindex,2));
        else
            CMap{rindex} = redblue(64);
        end
    end
    
    % Display Data
    axes(hA(rindex)); % set axes as current axes
    imagesc(data, CLim(rindex,:));
    colormap(hA(rindex), CMap{rindex});
    set(gca,'TickDir','out');
    hold on
    
    
    % Display spikes
    if showSpikes
        [Rindex, Tindex] = find(spikes>=spikeThresh);
        plot(Tindex, Rindex, 'g*');
    end
    
    % colorbar
    if showColorBar
        hcb = colorbar;
        ylabel(hcb, colorbarTitle);
        %         lab = get(hcb,'YTickLabel');
        %         set(hcb,'YTickLabel', cellstr(num2str(cellfun(@str2num,lab)+1)));
    end
    
    % Plot block lines
    if strcmp(sortType, 'stim')
        blockLines = cumsum(blockLength(1:end-1));
        plot(repmat([0,size(data,2)+1],length(blockLines),1)',repmat(blockLines+0.5,1,2)','k-','LineWidth',lineWidth); % trial dividing lines
    end
    
    % Plot stim lines
    plot(repmat([StimFrames(:,1)-.5, StimFrames(:,2)+.5]+xshift,2,1), repmat((0.5:1:size(StimFrames,1))',2,2), 'k--','LineWidth',lineWidth);
    
    % Label y-axis
    if strcmp(sortType, 'stim')
        ticks = cumsum(blockLength)-diff(cumsum([0;blockLength]))/2 + 0.5;
        set(hA(rindex),'YTick',ticks, 'YTickLabel',XLabels);
        if labelAxes
            ylabel('Stimulus');
        end
    elseif labelAxes
        ylabel('Trial')
    end
    
    % Label x-axis
    set(hA(rindex), 'XTick', 0.5:frameRate:framelimit(2)+0.5, 'XTickLabel', cellstr(num2str((0:framelimit/frameRate)')));
    if labelAxes
        xlabel(xlab);
    end
    
%     % Fix labels
%     if isempty(rois(ROIindex(rindex)).label)
%         rois(ROIindex(rindex)).label = {'unlabeled'};
%     end
%     
%     % Save plot to PDF
%     if saveToPDF
%         title(sprintf('ROI: %d, Label: %s', ROIindex(rindex), rois(ROIindex(rindex)).label{1}));
%         export_fig(hf, saveFile, '-append');
%         close(hf);
%     end
    
end %roi cycle

