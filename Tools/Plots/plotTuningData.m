function [h, hA, hF, normalize] = plotTuningData(Data, varargin)

ROIindex = [1,inf];

saveOut = false;
saveFile = ''; % filename to save plots to
StimIndex = [1,inf];

% Items to display
showDataPoints = false;
showStimStars = true;
showN = false;
showPValues = false;
ErrorBars = 'SE'; % 'SE', 'std', 'var', or '95CI'

% Plot colors & display options
normalize = false;
Colors = {'r'};
labelAxes = true;
XTick = [];
XTickLabel = {};
XTickLabelRotation = 90;
YLabel = 'Average Stimulus-Evoked dF/F';

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
if ~exist('Data', 'var') || isempty(Data)
    [Data, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Data)
        return
    end
    Data = fullfile(p, Data);
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1}
                index = index + 2;
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
            case 'ErrorBars'
                ErrorBars = varargin{index+1};
                index = index + 2;
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
            case 'YLabel'
                YLabel = varargin{index+1};
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
if ischar(Data)
    Data = {Data};
end
if iscellstr(Data)
    ROIFiles = Data;
    if saveOut && isempty(saveFile)
        [p,fn,~] = fileparts(ROIFiles{1});
        saveFile = fullfile(p, strcat(fn, '.pdf'));
    end
    Data = [];
    for findex = 1:numel(ROIFiles)
        load(ROIFiles{findex}, 'ROIdata', '-mat')
        Data = cat(1, Data, gatherROIdata(ROIdata,'Raw'));
    end
elseif isstruct(Data)
    Data = gatherROIdata(Data,'Raw');
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end


%% Determine ROIs to plot
if ischar(ROIindex) && strcmp(ROIindex, 'all')
    ROIindex = 1:size(Data,1);
elseif ROIindex(end) == inf
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:size(Data,1));
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
    StimIndex = [StimIndex(1:end-1),StimIndex(end-1)+1:size(Data,2)];
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
label=repmat({'none'},max(ROIindex),1);
h = nan(numStimuli,numel(ROIindex));
for index = 1:numel(ROIindex)
    rindex = ROIindex(index);
    
%     % Fix Label
%     if ~isfield(rois, 'label') || isempty(rois(rindex).label)
%         rois(rindex).label = {'none'};
%     end
    
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
            plot(s*ones(numel(Data{rindex,StimIndex(s)}),1), Data{rindex,StimIndex(s)}, 'k.') %plot raw data points for all stimuli
        end
    end
    
    % Gather data & compute error
    if iscell(Data)
        raw = Data(rindex,StimIndex);
        N = cellfun(@numel,raw);
        data = cellfun(@mean,raw);
    else
        data = Data(rindex,StimIndex);
        ErrorBars = 'none';
    end
    switch ErrorBars
        case 'SE'
            se = cellfun(@std,raw)./sqrt(N);
        case 'std'
            se = cellfun(@std,raw);
        case 'var'
            se = cellfun(@var,raw);
        case '95CI'
            se = 1.96*cellfun(@std,raw); % assumes normally distributed
    end
    
    % Normalize data
    if normalize(index)
        if isequal(normalize(index),1)
            normalize(index) = max(abs(data(2:end)));
        end
        data = data/normalize(index);
        se = se/normalize(index);
    end
    
    % Plot tuning
    for s = 1:numStimuli
        if ~strcmpi(ErrorBars,'none')
            h(s,index) = errorbar(s,data(s),se(s),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize); %plot curve
        else
            h(s,index) = plot(s,data(s),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize);
        end
    end
        
    % Set axes: limits, ticks, & labels
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
        ylabel(YLabel);
    end
    
    % Set Title
    if isequal(Title{AxesIndex(index)}, true)
        Title{AxesIndex(index)} = sprintf('ROI: %d, Label: %s',rindex,label{rindex});
    end
    if ~isempty(Title{AxesIndex(index)})
        title(Title{AxesIndex(index)});
    end
    
%     % Plot stars for stimuli that evoked a significant response
%     if showStimStars
%         Ydim = get(gca,'YLim');
%         for s = 1:numStimuli
%             if rois(rindex).PValueCorrected(StimIndex(s))<.05
%                 text(s,Ydim(2)-(Ydim(2)-Ydim(1))/10,'*','Color',[0,1,1],'FontSize',15,'HorizontalAlignment','center'); %display significance star
%             end
%         end
%     end
    
    % Display number of trials per average
    if showN
        Ydim = get(gca,'YLim');
        for s = 1:numStimuli
            text(s,Ydim(1)+(Ydim(2)-Ydim(1))/10,sprintf('n=%d',N(s)),'HorizontalAlignment','center');
        end
    end
    
%     % Display p-values
%     if showPValues
%         Ydim = get(gca,'YLim');
%         for s = 1:numStimuli
%             text(s,Ydim(1)+(Ydim(2)-Ydim(1))/20,sprintf('p=%.3f',rois(rindex).PValueCorrected(StimIndex(s))),'HorizontalAlignment','center');
%         end
%     end
    
    % Plot lines
    if ~isempty(HorzLines)
        XLim = get(gca,'XLim');
        for l = 1:numel(HorzLines)
            plot(XLim, [HorzLines(l) HorzLines(l)], 'b-');
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


