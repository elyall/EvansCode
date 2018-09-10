function [hA, CLim, CMap] = plotRaster(Data, DataInfo, varargin)

% Data to display
datatype = 'dFoF'; %'dFoF'; % 'dFoF' or 'raw'
sorttype = 'stim'; % 'stim' or 'trial' (will not show averages)
showdata = true;
showavgs = false;
TrialIndex = [1 inf];
StimOrder = [];

% Computation parameters
numSTDsOutlier = inf; % 'inf' means no trials are thrown out
ControlID = 0; % StimID of control stimulus ('nan' means no control trial)

% Display parameters
zscoreby = 'none'; %'all' or 'trial' or 'none'
avgPixelHeight = 6; % number of pixels of average plots height
frameRate = 15.45; %frame rate of the data
xlab='Times (s)'; %default x-axis label
framelimit = 'all'; % xlim (indexes data points)
lineWidth = .5;
ylabels = {};
showStimLines = false;
verticalLines = [];

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
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Trials','trials','TrialIndex'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'StimOrder'
                StimOrder = varargin{index+1};
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
            case 'showStimLines'
                showStimLines = true;
                index = index + 1;
            case 'lineWidth'
                lineWidth = varargin{index+1};
                index = index + 2;
            case 'verticalLines'
                verticalLines = varargin{index+1};
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
            case 'ControlID'
                ControlID = varargin{index+1};
                index = index + 2;
            case 'ylabels'
                ylabels = varargin{index+1};
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


%% Determine trials to plot
if TrialIndex(end)==inf
    TrialIndex = cat(2, TrialIndex(end-1), TrialIndex(end-1)+1:size(Data, 1));
end
if ~islogical(TrialIndex)
    TrialIndex = ismember(DataInfo.TrialIndex,TrialIndex);
end
numTrials = nnz(TrialIndex);


%% Determine stimulus of each trial
StimIndex = DataInfo.StimID;
if isempty(StimOrder)
    StimOrder = unique(StimIndex);
    temp = find(ismember(StimOrder,ControlID))'; % locate control stimuli
    StimOrder = StimOrder([temp, setdiff(1:numel(StimOrder),temp)]); % place control stimuli at front
end
if size(StimOrder,1)~=1
    StimOrder = StimOrder';
end
numStimuli = numel(StimOrder);
StimulusFrames = [repmat(DataInfo.numFramesBefore, numTrials, 1), DataInfo.numFramesBefore + DataInfo.numStimFrames(TrialIndex)];
StimulusFrames(ismember(DataInfo.StimID,ControlID),:) = nan; % do not show stim bars for control trials


%% Prepare data
numTrials = size(Data,1);

% Throw out unwanted trials
Data(~TrialIndex,:) = [];
StimIndex(~TrialIndex) = [];

% Throw out unwanted stimuli
Data(~ismember(StimIndex,StimOrder),:) = [];
StimIndex(~ismember(StimIndex,StimOrder)) = [];

% Sort stimuli
switch sorttype
    case 'stim'
        orderDict = bsxfun(@eq, StimIndex, StimOrder);
        [~,orderDict] = max(orderDict, [], 2);
        [~, displayIndex] = sort(orderDict);
        if isempty(ylabels)
            ylabels = cellstr(num2str(StimOrder'));
            [ylabels{ismember(StimOrder,ControlID)}] = deal('no contact');
        end
        blockLength = [];
        for sindex = 1:numStimuli
            blockLength = cat(1,blockLength,nnz(orderDict==sindex));
        end
    case 'trial'
        displayIndex = 1:size(Data,1);
        showavgs = false;
end

% Sort data
Data = Data(displayIndex, :);
StimIndex = StimIndex(displayIndex);
StimFrames = StimulusFrames(displayIndex, :);
if showSpikes
    spikes = spikes(displayIndex, :);
end

% Remove outliers for each stimulus
if strcmp(sorttype, 'stim') && (strcmp(zscoreby, 'all') || (strcmp(datatype, 'dFoF') && strcmp(zscoreby, 'none')))
    for sindex = 1:numStimuli
        n = inf;
        currentData = Data(StimIndex==sindex,:);
        [H,W] = size(currentData);
        currentData = reshape(zscore(currentData(:)), H, W);
        while n ~= blockLength(sindex)
            n = blockLength(sindex);
            bad = find(max(abs(currentData),[],2) > numSTDsOutlier); %throw out outlier trials
            if ~isempty(bad)
                currentData(bad, :) = [];
                blockLength(sindex) = blockLength(sindex) - numel(bad);
                bad = bad + sum(blockLength(1:sindex-1));
                Data(bad, :) = [];
                StimFrames(bad, :) = [];
                StimIndex(bad, :) = [];
                if showSpikes
                    spikes(bad, :) = [];
                end
            end
        end
    end
end

if showavgs
    % Compute average stimuli
    AvgStim = zeros(numStimuli*avgPixelHeight, size(Data, 2));
    AvgStimuliFrames = [];
    for sindex = 1:numStimuli
        AvgStim((sindex-1)*avgPixelHeight+1:(sindex-1)*avgPixelHeight+avgPixelHeight, :) = repmat(mean(Data(StimIndex==sindex, :), 1), avgPixelHeight, 1);
        AvgStimuliFrames = cat(1, AvgStimuliFrames, repmat(mode(StimFrames(StimIndex==sindex, :), 1), avgPixelHeight, 1));
    end
    
    % Add average stimuli to matrix to display
    if showdata
        Data = cat(1, Data, AvgStim);
        blockLength = cat(1, blockLength, repmat(avgPixelHeight, numStimuli,1));
        StimFrames = cat(1, StimFrames, AvgStimuliFrames);
        if ~isnan(ControlID)
            ylabels = cat(1, ylabels, [{'no contact avg'}; cellstr([num2str(StimIndex(~ismember(StimIndex, ControlID))), repmat(' avg', numStimuli-1,1)])]);
        else
            ylabels = cat(1, ylabels, cellstr([num2str(StimIndex), repmat(' avg', numStimuli,1)]));
        end
        %             labels = cat(1, labels, [{'control avg'}; cellstr([num2str((1:numStimuli-1)'),repmat(' avg',numStimuli-1,1)])]);
    else
        Data = AvgStim;
        blockLength = repmat(avgPixelHeight, numStimuli,1);
        StimFrames = AvgStimuliFrames;
        if ~isnan(ControlID)
            ylabels = [{'no contact'}, cellstr([num2str(setdiff(StimOrder,ControlID)'), repmat(' avg', numStimuli-1,1)])'];
        else
            ylabels = cellstr([num2str(StimIndex), repmat(' avg', numStimuli,1)]);
        end
        %             labels = [{'control avg'}; cellstr([num2str((1:numStimuli-1)'),repmat(' avg',numStimuli-1,1)])];
    end
end

% Z-score data
switch zscoreby
    case 'trial'
        Data = zscore(Data,0,2); %zscore each row individually
        colorbarTitle = 'z-score (per trial)';
    case 'all'
        [H,W] = size(Data);
        Data = reshape(zscore(Data(:)),H,W);
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
    Data = Data(:,framelimit(1):framelimit(2));
    if showSpikes
        spikes = spikes(:,framelimit(1):framelimit(2));
    end
else
    xshift = 0;
end


%% Create axes
if isempty(hA)
    figure(); hA = axes();
else
    axes(hA); % set axes as current axes
end

% Determine color attritbutes
if isempty(CLim)
    CLim = [nanmin(Data(:)), nanmax(Data(:))];
end
if isempty(CMap)
%     if ~strcmp(zscoreby, 'none') || strcmp(datatype, 'dFoF')
        CMap = b2r(CLim(1),CLim(2));
%     else
%         CMap = redblue(256);
%         CMap = parula(256);
%         CMap = hot(256);
%     end
end


%% Display Data
imagesc(Data, CLim);
colormap(hA, CMap);
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
    plot(repmat([0,size(Data,2)+1],length(blockLines),1)',repmat(blockLines+0.5,1,2)','k-','LineWidth',lineWidth); % trial dividing lines
end

% Plot stim lines
if showStimLines
    plot(repmat(StimFrames(:,1)-.5+xshift,2,1), repmat((0.5:1:size(StimFrames,1))',2,1), 'k--','LineWidth',lineWidth);
    plot(repmat([StimFrames(:,1)-.5, StimFrames(:,2)+.5]+xshift,2,1), repmat((0.5:1:size(StimFrames,1))',2,2), 'k--','LineWidth',lineWidth);
end

% Plot user-specified vertical lines
if ~isempty(verticalLines)
    YLim = get(gca,'YLim');
    for index = 1:numel(verticalLines)
        plot([verticalLines(index),verticalLines(index)],YLim, 'k--','LineWidth',lineWidth);
    end
end

% Label y-axis
if strcmp(sorttype, 'stim')
    ticks = cumsum(blockLength)-diff(cumsum([0;blockLength]))/2 + 0.5;
    set(hA,'YTick',ticks, 'YTickLabel',ylabels);
    ylabel('Stimulus');
else
    ylabel('Trial')
end

% Label x-axis
set(hA, 'XTick', 0.5:frameRate:framelimit(2)+0.5, 'XTickLabel', cellstr(num2str((0:framelimit/frameRate)')));
xlabel(xlab);


