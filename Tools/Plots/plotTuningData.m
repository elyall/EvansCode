function [h, hA, hF, normalize] = plotTuningData(Data, varargin)

ROIindex = [1,inf];

saveOut = false;
saveFile = ''; % filename to save plots to

StimIndex = [1,inf];
Grouping = []; % cell array that where each cell contains x-coordinates of data points corresponding to that line object 
PValues = [];

% Items to display
showDataPoints = false;
showStimStars = false;
showN = false;
showPValues = false;
ErrorBars = '95CI'; % 'SE', 'std', 'var', '95CI', 'none', or matrix of size [1or2,numStimuli,nROIs]

% Plot colors & display options
PlotType = 'normal'; % 'normal' or 'polar'
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
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case 'Grouping'
                Grouping = varargin{index+1};
                index = index + 2;
            case 'PValues'
                PValues = varargin{index+1};
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
            case 'PlotType'
                PlotType = varargin{index+1};
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
            case {'XTickRot','XTickLabelRotation'}
                XTickLabelRotation = varargin{index+1};
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

if isempty(PValues)
    showStimStars = false;
    showPValues = false;
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

if iscell(Data) && showDataPoints && any(cellfun(@isrow,Data)) % needs to be column for cat() call
    Data = cellfun(@transpose,Data,'UniformOutput',false);
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
if isempty(Grouping)
    Grouping = mat2cell(1:numStimuli,1,ones(1,numStimuli));
end
if strcmp(PlotType,'polar') % set last data point to 360
    theta = linspace(0,360,numStimuli+1);
    theta = theta(1:end-1);
end

%% Determine plotting colors
if ischar(Colors)
    Colors = repmat({Colors}, numStimuli, 1);
elseif isnumeric(Colors)
    if size(Colors,1)==1
        Colors = repmat({Colors}, numStimuli, 1);
    else
        Colors = mat2cell(Colors,ones(size(Colors,1),1),size(Colors,2));
    end
elseif iscell(Colors) && numel(Colors) == 1
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
        switch PlotType
            case 'normal'
                hA(AxesIndex(index)) = axes();
            case 'polar'
                hA(AxesIndex(index)) = polaraxes();
        end
    else
        axes(hA(AxesIndex(index)));
        hF = get(hA(AxesIndex(index)), 'Parent');
    end
    hold on
    
    % Gather data
    if iscell(Data)
        raw = Data(rindex,StimIndex);
        N = cellfun(@numel,raw);
        data = cellfun(@mean,raw);
    else
        data = Data(rindex,StimIndex);
    end
    
    % Determine error
    if ischar(ErrorBars)
        switch ErrorBars
            case 'none'
                se = [];
            case 'SE'
                se = cellfun(@std,raw)./sqrt(N);
            case 'std'
                se = cellfun(@std,raw);
            case 'var'
                se = cellfun(@var,raw);
            case '95CI'
                se = 1.96*cellfun(@std,raw)./sqrt(N); % assumes normally distributed
        end
    elseif isnumeric(ErrorBars)
        se = squeeze(ErrorBars(:,StimIndex,rindex));
    else
        ErrorBars = 'none';
    end
    if ~isempty(PValues)
        p = PValues(rindex,StimIndex);
    end
    
    % Normalize data
    if normalize(index)
        if isequal(normalize(index),1)
            normalize(index) = max(abs(data(2:end)));
        end
        data = data/normalize(index);
        se = se/normalize(index);
    end
    
    switch PlotType
        
        case 'normal'
            
            % Plot tuning
            for s = 1:numel(Grouping)
                if isempty(se)
                    h(s,index) = plot(Grouping{s},data(Grouping{s}),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize);
                elseif size(se,1)==1
                    h(s,index) = errorbar(Grouping{s},data(Grouping{s}),se(Grouping{s}),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize);
                else
                    h(s,index) = errorbar(Grouping{s},data(Grouping{s}),se(1,Grouping{s}),se(2,Grouping{s}),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize);
                end
            end
            
            % Plot each trial's average dF/F for each stimulus
            if showDataPoints
                plotSpread(cat(1,raw{:}),'distributionIdx',repelem(1:numStimuli,cellfun(@numel,raw))','binWidth',0.4,'distributionColors',[.5,.5,.5]) %plot raw data points for all stimuli        for s = 1:numStimuli
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
            
            % Plot stars for stimuli that evoked a significant response
            if showStimStars
                Ydim = get(gca,'YLim');
                for s = find(p<.05)
                    text(s,Ydim(2)-(Ydim(2)-Ydim(1))/10,'*','Color',[0,1,1],'FontSize',15,'HorizontalAlignment','center'); %display significance star
                end
            end
            
            % Display number of trials per average
            if showN
                Ydim = get(gca,'YLim');
                for s = 1:numStimuli
                    text(s,Ydim(1)+(Ydim(2)-Ydim(1))/10,sprintf('n=%d',N(s)),'HorizontalAlignment','center');
                end
            end
            
            % Display p-values
            if showPValues
                Ydim = get(gca,'YLim');
                for s = 1:numStimuli
                    text(s,Ydim(1)+(Ydim(2)-Ydim(1))/20,sprintf('p=%.3f',p(s)),'HorizontalAlignment','center');
                end
            end
            
            % Plot lines
            if ~isempty(HorzLines)
                XLim = get(gca,'XLim');
                for l = 1:numel(HorzLines)
                    plot(XLim, [HorzLines(l) HorzLines(l)], 'k-');
                end
                temp = get(gca,'Children');
                set(gca,'Children',temp([2:end,1]));
                set(gca,'XLim',XLim);
            end
            if ~isempty(VertLines)
                YLim = get(gca,'YLim');
                for l = 1:numel(VertLines)
                    plot([VertLines(l) VertLines(l)], YLim, 'k--');
                end
                set(gca,'YLim',YLim);
            end
            
        case {'polar','circ','circular'}
            
            % Plot tuning
            for s = 1:numel(Grouping)
%                 if isempty(se)
                    h(s,index) = polarplot(deg2rad(theta(Grouping{s})),data(Grouping{s}),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize);
%                 elseif size(se,1)==1
%                     h(s,index) = errorbar(Grouping{s},data(Grouping{s}),se(Grouping{s}),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize);
%                 else
%                     h(s,index) = errorbar(Grouping{s},data(Grouping{s}),se(1,Grouping{s}),se(2,Grouping{s}),'Color',Colors{s},'LineStyle','-','LineWidth',BarWidth,'Marker','.','MarkerSize',MarkerSize);
%                 end
            end
            
            % Plot each trial's average dF/F for each stimulus
            if showDataPoints
                for s = 1:numel(raw)
                    polarscatter(deg2rad(theta(s))*ones(numel(raw{s}),1),raw{s},5,[.5,.5,.5],'.');
                end
            end
            
            % Set axes: limits, ticks, & labels
            if ~isempty(YLim)
                if isnumeric(YLim)
                    rlim(YLim);
                elseif ischar(YLim)
                    axis(YLim);
                end
            end
            if isempty(XTick)
                thetaticks(theta);
            else
                thetaticks(XTick);
            end
            if ~isempty(XTickLabel)
                thetaticklabels(XTickLabel);
            end
            % thetalim([0,360]);
            % set(gca,'XTickLabelRotation',XTickLabelRotation);
            % if labelAxes
            %     xlabel('Stimulus');
            %     ylabel(YLabel);
            % end
            
%             % Plot stars for stimuli that evoked a significant response
%             if showStimStars
%                 Ydim = get(gca,'YLim');
%                 for s = find(p<.05)
%                     text(s,Ydim(2)-(Ydim(2)-Ydim(1))/10,'*','Color',[0,1,1],'FontSize',15,'HorizontalAlignment','center'); %display significance star
%                 end
%             end
%             
%             % Display number of trials per average
%             if showN
%                 Ydim = get(gca,'YLim');
%                 for s = 1:numStimuli
%                     text(s,Ydim(1)+(Ydim(2)-Ydim(1))/10,sprintf('n=%d',N(s)),'HorizontalAlignment','center');
%                 end
%             end
%             
%             % Display p-values
%             if showPValues
%                 Ydim = get(gca,'YLim');
%                 for s = 1:numStimuli
%                     text(s,Ydim(1)+(Ydim(2)-Ydim(1))/20,sprintf('p=%.3f',p(s)),'HorizontalAlignment','center');
%                 end
%             end
%             
%             % Plot lines
%             if ~isempty(HorzLines)
%                 XLim = get(gca,'XLim');
%                 for l = 1:numel(HorzLines)
%                     plot(XLim, [HorzLines(l) HorzLines(l)], 'k-');
%                 end
%                 temp = get(gca,'Children');
%                 set(gca,'Children',temp([2:end,1]));
%                 set(gca,'XLim',XLim);
%             end
%             if ~isempty(VertLines)
%                 YLim = get(gca,'YLim');
%                 for l = 1:numel(VertLines)
%                     plot([VertLines(l) VertLines(l)], YLim, 'k--');
%                 end
%                 set(gca,'YLim',YLim);
%             end
            
    end %switch PlotType
    
    % Set Title
    if isequal(Title{AxesIndex(index)}, true)
        Title{AxesIndex(index)} = sprintf('ROI: %d, Label: %s',rindex,label{rindex});
    end
    if ~isempty(Title{AxesIndex(index)})
        title(Title{AxesIndex(index)});
    end
            
    hold off
    
    % Save plot to PDF
    if saveOut && ismember(index, lastROI)
        drawnow;
        export_fig(hF, saveFile, '-append');
        close(hF);
    end
    
end %cycle ROIs


