function [hA, CLim, CMap] = plotIndividualRaster(rois, ROIindex, DataInfo, varargin)

saveToPDF = false;
saveFile = 'test.pdf'; %filename to save plots to

% Data to display
datatype = 'dFoF'; %'dFoF'; % 'dFoF' or 'raw'
sorttype = 'stim'; % 'stim' or 'trial' (will not show averages)
showdata = true;
showavgs = true;
TrialIndex = [1 inf];
StimIndex = [];

% Computation parameters
minrunspeed = 0; % (will only run if 'RunningSpeed' found in ROIFile)
numSTDsOutlier = inf; % Inf means no trials are thrown out

% Display parameters
zscoreby = 'none'; %'all' or 'trial' or 'none'
avgPixelHeight = 6; % number of pixels of average plots height
frameRate = 15.45; %frame rate of the data
xlab='Frames'; %default x-axis label
framelimit = 'all'; % xlim (indexes data points)
lineWidth = .5;

% Color properties
showColorBar = false;
CLim = [];
CMap = {};

% Spike parameters
showSpikes = false;
spikeThresh = .18;

% Constant parameters
controlindex = nan; %index of control stimulus ID relative to other StimIDs (ID is minimum value => 1)
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
    ROIindex = 'all';
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Trials','trials','TrialIndex'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case 'datatype'
                datatype = varargin{index+1};
                index = index + 2;
            case 'sorttype'
                sorttype = varargin{index+1};
                index = index + 2;
            case 'showdata'
                showdata = varargin{index+1};
                index = index + 2;
            case 'showavgs'
                showavgs = varargin{index+1};
                index = index + 2;
            case 'minrunspeed'
                minrunspeed = varargin{index+1};
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
            case 'saveToPDF'
                saveToPDF = true;
                index = index + 1;
            case 'filename'
                saveFile = varargin{index+1};
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


if ~exist('RunningSpeed', 'var') && minrunspeed ~= 0
    minrunspeed = 0;
    warning('RunningSpeed info does not exist in ROIFile. Will not remove any non-running trials.');
end


%% Determine trials to plot
if TrialIndex(end)==inf
    TrialIndex = cat(2, TrialIndex(end-1), TrialIndex(end-1)+1:size(rois(1).data, 1));
elseif islogical(TrialIndex)
    TrialIndex = find(TrialIndex);
end
numTrials = numel(TrialIndex);


%% Determine stimulus of each trial
if isempty(StimIndex)
    [StimIDs, ~, TrialID] = unique(DataInfo.StimID(TrialIndex));
else
    if numel(StimIndex) == numTrials
        [StimIDs, ~, TrialID] = unique(StimIndex);
    else
        [StimIDs, ~, TrialID] = unique(StimIndex(TrialIndex));
    end
end
numStimuli = numel(StimIDs);
StimulusFrames = [repmat(DataInfo.numFramesBefore, numTrials, 1), DataInfo.numFramesBefore + DataInfo.numStimFrames(TrialIndex)];
StimulusFrames(TrialID==controlindex,:) = nan; % do not show stim bars for control trials


%% Determine ROIs to plot
if ischar(ROIindex) && strcmp(ROIindex, 'all')
    ROIindex = 1:numel(rois);
elseif ROIindex(end) == inf
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
    switch datatype
        case 'dFoF'
            data = rois(ROIindex(rindex)).dFoF(TrialIndex, :);
        case 'raw'
            data = rois(ROIindex(rindex)).data(TrialIndex, :);
    end
    if showSpikes
        spikes = rois(ROIindex(rindex)).spikes(TrialIndex,:);
    end
    
    % Remove non-running trials
    if minrunspeed
        data(NotRunning, :) = [];
        if showSpikes
            spikes(NotRunning, :) = [];
        end
    end
    
    % Sort stimuli
    switch sorttype
        case 'stim'
            TrialID(TrialID == controlindex) = -Inf; %this doesn't matter if controlindex is the smallest value (but does if controlindex is higher)
            [TrialID, displayIndex] = sort(TrialID);
            TrialID(TrialID == -Inf) = controlindex;
            labels = cellstr(num2str(StimIDs));
%             labels = [{'control'}; cellstr(num2str((1:numStimuli-1)'))];
            blockLength = [];
            for sindex = 1:numStimuli
                blockLength = cat(1,blockLength,sum(TrialID==sindex));
            end
        case 'trial'
            displayIndex = 1:numTrials;
            showavgs = false;
    end
    
    % Sort data
    data = data(displayIndex, :);
    StimFrames = StimulusFrames(displayIndex, :);
    if showSpikes
        spikes = spikes(displayIndex, :);
    end
    
    % Remove outliers for each stimulus
    if strcmp(sorttype, 'stim') && (strcmp(zscoreby, 'all') || (strcmp(datatype, 'dFoF') && strcmp(zscoreby, 'none')))
        for sindex = 1:numStimuli
            n = inf;
            currentData = data(TrialID==sindex,:);
            [H,W] = size(currentData);
            currentData = reshape(zscore(currentData(:)), H, W);
            while n ~= blockLength(sindex)
                n = blockLength(sindex);
                bad = find(max(abs(currentData),[],2) > numSTDsOutlier); %throw out outlier trials
                if ~isempty(bad)
                    currentData(bad, :) = [];
                    blockLength(sindex) = blockLength(sindex) - numel(bad);
                    bad = bad + sum(blockLength(1:sindex-1));
                    data(bad, :) = [];
                    StimFrames(bad, :) = [];
                    TrialID(bad, :) = [];
                    if showSpikes
                        spikes(bad, :) = [];
                    end
                end
            end
        end
    end
        
    if showavgs
        % Compute average stimuli
        AvgStim = zeros(numStimuli*avgPixelHeight, size(data, 2));
        AvgStimuliFrames = [];
        for sindex = 1:numStimuli
            AvgStim((sindex-1)*avgPixelHeight+1:(sindex-1)*avgPixelHeight+avgPixelHeight, :) = repmat(mean(data(TrialID==sindex, :), 1), avgPixelHeight, 1);
            AvgStimuliFrames = cat(1, AvgStimuliFrames, repmat(mode(StimFrames(TrialID==sindex, :), 1), avgPixelHeight, 1));
        end
        
        % Add average stimuli to matrix to display
        if showdata
            data = cat(1, data, AvgStim);
            blockLength = cat(1, blockLength, repmat(avgPixelHeight, numStimuli,1));
            StimFrames = cat(1, StimFrames, AvgStimuliFrames);
            labels = cat(1, labels, cellstr([num2str(StimIDs), repmat(' avg', numStimuli,1)]));
%             labels = cat(1, labels, [{'control avg'}; cellstr([num2str((1:numStimuli-1)'),repmat(' avg',numStimuli-1,1)])]);
        else
            data = AvgStim;
            blockLength = repmat(avgPixelHeight, numStimuli,1);
            StimFrames = AvgStimuliFrames;
            labels = cellstr([num2str(StimIDs), repmat(' avg', numStimuli,1)]);
%             labels = [{'control avg'}; cellstr([num2str((1:numStimuli-1)'),repmat(' avg',numStimuli-1,1)])];
        end
    end
    
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
            switch datatype
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
    if ~strcmp(zscoreby, 'none') || strcmp(datatype, 'dFoF')
        CMap{rindex} = b2r(CLim(rindex,1),CLim(rindex,2));
    else
        CMap{rindex} = redblue(64);
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
    if strcmp(sorttype, 'stim')
        blockLines = cumsum(blockLength);
        blockLines = blockLines(1:end-1);
        plot(repmat([0,size(data,2)+1],length(blockLines),1)',repmat(blockLines+0.5,1,2)','k-','LineWidth',lineWidth); % trial dividing lines
    end
    
    % Plot stim lines
    plot(repmat(StimFrames(:,1)-.5+xshift,2,1), repmat((0.5:1:size(StimFrames,1))',2,1), 'k--','LineWidth',lineWidth);
%     plot(repmat([StimFrames(:,1)-.5, StimFrames(:,2)+.5]+xshift,2,1), repmat((0.5:1:size(StimFrames,1))',2,2), 'k--','LineWidth',lineWidth);
    
    % Label y-axis
    if strcmp(sorttype, 'stim')
        ticks = cumsum(blockLength)-diff(cumsum([0;blockLength]))/2 + 0.5;
        set(hA(rindex),'YTick',ticks, 'YTickLabel',labels);
        ylabel('Stimulus');
    else
        ylabel('Trial')
    end
    
    % Label x-axis
    set(hA(rindex), 'XTick', 0.5:frameRate:framelimit(2)+0.5, 'XTickLabel', cellstr(num2str((0:framelimit/frameRate)')));
    xlabel(xlab);
    
    % Fix labels
    if isempty(rois(ROIindex(rindex)).label)
        rois(ROIindex(rindex)).label = {'unlabeled'};
    end
    
    % Save plot to PDF
    if saveToPDF
        title(sprintf('ROI: %d, Label: %s', ROIindex(rindex), rois(ROIindex(rindex)).label{1}));
        export_fig(hf, saveFile, '-append');
        close(hf);
    end
    
end %roi cycle

