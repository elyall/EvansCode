function [h, hA, hF, normalize] = plotTuning(rois, ROIindex, varargin)

saveOut = false;
saveFile = ''; % filename to save plots to
StimIndex = [1,inf];

% Items to display
showDataPoints = false;
showStimStars = true;
showN = false;
showPValues = false;

% Plot colors & display options
normalize = false;
Colors = {'r'};
labelAxes = true;
XTick = [];
XTickLabel = {};
XTickLabelRotation = 90;

% Lines
HorzLines = [];
VertLines = [];
LineStyle = '--';
LineColor = 'k';
LineWidth = 2;

% Constant variables
AxesIndex = [];
hA = [];
Title = true;
YLim = [];
MarkerSize = 20; % 20
BarWidth = 2;

directory = cd;

%% Parse input arguments
if ~exist('rois', 'var') || isempty(rois)
    [rois, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(rois)
        return
    end
    rois = fullfile(p, rois);
end

if ~exist('ROIindex','var') || isempty(ROIindex)
    ROIindex = 'all';
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'showDataPoints'
                showDataPoints = true;
                index = index + 1;
            case 'showStimStars'
                showStimStars = true;
                index = index + 1;
            case 'showN'
                showN = true;
                index = index + 1;
            case 'showPValues'
                showPValues = true;
                index = index + 1;
            case 'normalize'
                normalize = varargin{index+1};
                index = index + 2;
            case 'Colors'
                Colors = varargin{index+1};
                index = index + 2;
            case 'labelAxes'
                labelAxes = true;
                index = index + 1;
            case 'XTick'
                XTick = varargin{index+1};
                index = index + 2;
            case 'XTickLabel'
                XTickLabel = varargin{index+1};
                index = index + 2;
            case 'HorzLines'
                HorzLines = varargin{index+1};
                index = index + 2;
            case 'VertLines'
                VertLines = varargin{index+1};
                index = index + 2;
            case 'LineStyle'
                LineStyle = varargin{index+1};
                index = index + 2;
            case 'LineColor'
                LineColor = varargin{index+1};
                index = index + 2;
            case 'LineWidth'
                LineWidth = varargin{index+1};
                index = index + 2;
            case 'AxesIndex'
                AxesIndex = varargin{index+1};
                index = index + 2;
            case 'axes'
                hA = varargin{index+1};
                index = index + 2;
            case 'Title'
                Title = varargin{index+1};
                index = index + 2;
            case 'YLim'
                YLim = varargin{index+1};
                index = index + 2;
            case 'showZero'
                showZero = true;
                index = index + 1;
            case 'BarWidth'
                BarWidth = varargin{index+1};
                index = index + 2;
            case 'MarkerSize'
                MarkerSize = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
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
    if saveOut && isempty(saveFile)
        [p,fn,~] = fileparts(ROIFile);
        saveFile = fullfile(p, strcat(fn, '.pdf'));
    end
    load(ROIFile, 'ROIdata', '-mat');
    rois = ROIdata.rois;
elseif iscellstr(rois)
    ROIFiles = rois;
    if saveOut && isempty(saveFile)
        [p,fn,~] = fileparts(ROIFiles{1});
        saveFile = fullfile(p, strcat(fn, '.pdf'));
    end
    rois = [];
    for findex = 1:numel(ROIFiles)
        load(ROIFiles{findex}, 'ROIdata', '-mat')
        rois = cat(1, rois, ROIdata.rois);
    end
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end


%% Determine ROIs to plot
if ischar(ROIindex) && strcmp(ROIindex, 'all')
    ROIindex = 1:numel(rois);
elseif ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:numel(rois));
end
numROIs = numel(ROIindex);

if numel(normalize) == 1
    normalize = repmat(normalize,numROIs,1);
end


%% Determine number of axes
if ~isempty(hA) && isempty(AxesIndex)
    if numel(hA)==1
        AxesIndex = ones(1,numROIs); % plots all to same plot
    elseif numel(hA)==numROIs
        AxesIndex = 1:numROIs;
    end
elseif isempty(AxesIndex)
    AxesIndex = 1:numROIs;
    hA = nan(numROIs, 1);
elseif isempty(hA)
    [~,~,AxesIndex] = unique(AxesIndex); % in case AxesIndex is not contiguous
    hA = nan(max(AxesIndex), 1);
end
numAxes = numel(hA);

% Determine which ROI is last for each figure
[~, lastROI] = unique(AxesIndex, 'last');

% Initialize titles container
if isempty(Title)
    Title = cell(numAxes, 1);
elseif ~iscell(Title)
    Title = repmat({Title}, numAxes, 1);
end


%% Determine stimuli to display
if StimIndex(end) == inf
    StimIndex = [StimIndex(1:end-1),StimIndex(end-1)+1:numel(rois(1).curve)];
end
numStimuli = numel(StimIndex);


%% Determine plotting colors
if ischar(Colors)
    Colors = repmat({Colors}, numStimuli, 1);
elseif isnumeric(Colors)
    Colors = repmat({Colors}, numStimuli, 1);
elseif iscellstr(Colors) && numel(Colors) == 1
    Colors = repmat(Colors, numStimuli, 1);
end


%% Plot each specified ROI's tuning
h = nan(numStimuli,numel(ROIindex));
for index = 1:numel(ROIindex)
    rindex = ROIindex(index);
    
    % Fix Label
    if ~isfield(rois, 'label') || isempty(rois(rindex).label)
        rois(rindex).label = {'none'};
    end
    
    % Select axes
    if isnumeric(hA(AxesIndex(index))) && isnan(hA(AxesIndex(index)))
        hF = figure();
        hA(AxesIndex(index)) = axes();
    else
        axes(hA(AxesIndex(index)));
        hF = get(hA(AxesIndex(index)), 'Parent');
    end
    hold on
    
    % Plot each trial's average dF/F for each stimulus
    if showDataPoints
        for s = 1:numStimuli
            plot(s*ones(rois(rindex).nTrials(s),1), rois(rindex).Raw{StimIndex(s)}, 'k.') %plot raw data points for all stimuli
        end
    end
    
    % Plot tuning curves
    if normalize(index)
        if isequal(normalize(index),1)
            normalize(index) = max(abs(rois(rindex).curve(2:end)));
        end
        data = rois(rindex).curve/normalize(index);
        se = rois(rindex).StdError/normalize(index);
    else
        data = rois(rindex).curve;
        se = rois(rindex).StdError;
    end
    for s = 1:numStimuli
        h(s,index) = errorbar(s,data(StimIndex(s)),se(StimIndex(s)),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize); %plot curve
    end
        
    % Set axes labels
    if ~isempty(YLim)
        if isnumeric(YLim)
            ylim(YLim);
        elseif ischar(YLim)
            axis(YLim);
        end
    end
    if ~isempty(XTick) && isempty(XTickLabel)
        set(gca,'XTick',XTick);
    else
        set(gca,'XTick',XTick,'XTickLabel',XTickLabel)
    end
    xlim([0,numStimuli+1]);
    set(gca,'XTickLabelRotation',XTickLabelRotation);
    if labelAxes
        xlabel('Stimulus');
        ylabel('Average Stimulus-Evoked dF/F');
    end
    
    % Set Title
    if isequal(Title{AxesIndex(index)}, true)
        Title{AxesIndex(index)} = sprintf('ROI: %d, Label: %s',rindex,rois(rindex).label{1});
    end
    if ~isempty(Title{AxesIndex(index)})
        title(Title{AxesIndex(index)});
    end
    
    % Plot stars for stimuli that evoked a significant response
    if showStimStars
        Ydim = get(gca,'YLim');
        for s = 1:numStimuli
            if rois(rindex).PValueCorrected(StimIndex(s))<.05
                text(s,Ydim(2)-(Ydim(2)-Ydim(1))/10,'*','Color',[0,1,1],'FontSize',15,'HorizontalAlignment','center'); %display significance star
            end
        end
    end
    
    % Display number of trials per average
    if showN
        Ydim = get(gca,'YLim');
        for s = 1:numStimuli
            text(s,Ydim(1)+(Ydim(2)-Ydim(1))/10,sprintf('n=%d',rois(rindex).nTrials(StimIndex(s))),'HorizontalAlignment','center');
        end
    end
    
    % Display p-values
    if showPValues
        Ydim = get(gca,'YLim');
        for s = 1:numStimuli
            text(s,Ydim(1)+(Ydim(2)-Ydim(1))/20,sprintf('p=%.3f',rois(rindex).PValueCorrected(StimIndex(s))),'HorizontalAlignment','center');
        end
    end
    
    % Plot lines
    if ~isempty(HorzLines)
        XLim = get(gca,'XLim');
        for l = 1:numel(HorzLines)
            plot(XLim, [HorzLines(l) HorzLines(l)], 'k--');
        end
        set(gca,'XLim',XLim);
    end
    if ~isempty(VertLines)
        YLim = get(gca,'YLim');
        for l = 1:numel(VertLines)
            plot([VertLines(l) VertLines(l)], YLim, 'k--');
        end
        set(gca,'YLim',YLim);
    end
    
    hold off
    
    % Save plot to PDF
    if saveOut && ismember(index, lastROI)
        drawnow;
        export_fig(hF, saveFile, '-append');
        close(hF);
    end
    
end %cycle ROIs


