function [hA, hF, normalize] = plotTuningCurve(rois, ROIindex, varargin)

saveOut = false;
saveFile = ''; % filename to save plots to

% Items to display
showControl = true;
showFit = false;
showDataPoints = false;
showStimStars = false;
showN = false;
showPValues = false;
PWCZ = [];

% Plot colors & display options
normalize = false;
curveColor = 'r';
fitColor = 'g';
fitLegendLocation = 'NorthWest'; % 'NorthWest' or 'NorthEast'
fontSize = 12;
labelAxes = false;

% Constant variables
AxesIndex = [];
hA = [];
Title = true;
YLim = [];
Legend = {};
LegendLocation = 'NorthEast';
MarkerSize = 1; % 20
LineWidth = 1; % 2
showZero = false;

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
            case 'showControl'
                showControl = true;
                index = index + 1;
            case 'showFit'
                showFit = true;
                index = index + 1;
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
            case 'PWCZ'
                PWCZ = varargin{index+1};
                index = index + 2;
            case 'normalize'
                normalize = varargin{index+1};
                index = index + 2;
            case 'curveColor'
                curveColor = varargin{index+1};
                index = index + 2;
            case 'fitColor'
                fitColor = varargin{index+1};
                index = index + 2;
            case 'legendLocation'
                fitLegendLocation = varargin{index+1};
                index = index + 2;
            case 'fontSize'
                fontSize = varargin{index+1};
                index = index + 2;
            case 'labelAxes'
                labelAxes = true;
                index = index + 1;
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
            case 'Legend'
                Legend = varargin{index+1};
                index = index + 2;
            case 'LegendLocation'
                LegendLocation = varargin{index+1};
                index = index + 2;
            case 'LineWidth'
                LineWidth = varargin{index+1};
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
elseif ROIindex(end) == inf;
    ROIindex = cat(2, ROIindex(1:end-1), ROIindex(1:end-1)+1:numel(rois));
end
numROIs = numel(ROIindex);

if numel(normalize) == 1
    normalize = repmat(normalize,numROIs,1);
end


%% Determine number of axes
if isempty(AxesIndex)
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


%% Determine plotting colors
if ischar(curveColor)
    curveColor = repmat({curveColor}, numROIs, 1);
elseif isnumeric(curveColor)
    curveColor = repmat({curveColor}, numROIs, 1);
elseif iscellstr(curveColor) && numel(curveColor) == 1
    curveColor = repmat(curveColor, numROIs, 1);
end

if ischar(fitColor)
    fitColor = repmat({fitColor}, numROIs, 1);
elseif iscellstr(fitColor) && numel(fitColor) == 1
    fitColor = repmat(fitColor, numROIs, 1);
end


%% Plot each specified ROI's tuning

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
    numStimuli = numel(rois(rindex).curve);
    if showDataPoints
        if showControl
            plot(ones(rois(rindex).nTrials(1),1), rois(rindex).Raw{1}, 'k.') %plot raw data points for control position
        end
        for s = 2:numStimuli
            plot((s-~showControl)*ones(rois(rindex).nTrials(s),1), rois(rindex).Raw{s}, 'k.') %plot raw data points for all stimuli
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
    if showControl
        errorbar(1,data(1),se(1),'Color',curveColor{index},'LineStyle','-','LineWidth',LineWidth,'Marker','.','MarkerSize',MarkerSize);  %plot control position
    end
    errorbar(1+showControl:numStimuli-1+showControl,data(2:end),se(2:end),'Color',curveColor{index},'LineStyle','-','LineWidth',LineWidth,'Marker','.','MarkerSize',MarkerSize); %plot curve
    
    % Plot fit
    if showFit && numStimuli > 1
        yfit = feval(rois(rindex).Fit,1:0.001:numStimuli-1);
        plot(2:0.001:numStimuli, yfit + rois(rindex).offset, 'Color', fitColor{index}, 'LineStyle', '-','LineWidth',1.25);
        Ydim = get(gca,'YLim');
        Xdim = get(gca,'XLim');
        switch fitLegendLocation
            case 'NorthEast'
                text(Xdim(2),Ydim(2),sprintf('r^2=%.2f\nFWHM=%.2f', rois(rindex).rsquare, rois(rindex).Coeff(3)),'Color',fitColor{index},'FontSize',fontSize,'HorizontalAlignment','right','VerticalAlignment','top');
            case 'NorthWest'
                text(Xdim(1),Ydim(2),sprintf('r^2=%.2f\nFWHM=%.2f', rois(rindex).rsquare, rois(rindex).Coeff(3)),'Color',fitColor{index},'FontSize',fontSize,'HorizontalAlignment','left','VerticalAlignment','top');
        end
    end
    
    % Set axes labels
    if ~isempty(YLim)
        if isnumeric(YLim)
            ylim(YLim);
        elseif ischar(YLim)
            axis(YLim);
        end
    end
    if showControl
        set(gca,'XTick',1:numStimuli,'XTickLabel',[{'no contact'}; cellstr(num2str((1:numStimuli-1)'))]);
        xlim([0,numStimuli+1]);
    else
        set(gca,'XTick',1:numStimuli-1,'XTickLabel',[cellstr(num2str((1:numStimuli-1)'))]);
        xlim([0,numStimuli]);
    end
    if labelAxes
        xlabel('Position');
        ylabel('Average Stimulus-Evoked dF/F');
    end
    
    % Place legend
    if ~isempty(Legend) && ismember(index, lastROI)
        legend(Legend, 'Location', LegendLocation);
        legend boxoff
    end
    
    % Show PWCZ
    if ~isempty(PWCZ)
        temp = ylim(gca);
        patch([repmat(PWCZ(1)+showControl-0.5,2,1);repmat(PWCZ(end)+showControl+0.5,2,1)], [temp';flip(temp')],[.5,.5,.5],...
            'FaceAlpha',    0.2,...
            'EdgeColor',    [.5,.5,.5],...
            'EdgeAlpha',    0.2);
        temp = get(gca,'children');             % get ordering
        set(gca,'children',temp([2:end,1]));    % move patch to back of plot 
        % plot(repmat([PWCZ(1)+showControl-0.5, PWCZ(end)+showControl+0.5],2,1), repmat(temp',1,2), 'k--'); % line boundaries
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
        sigindices = find(rois(rindex).SigPos); %locate any responses that are significant
        if ~isempty(sigindices)
            text(sigindices,(Ydim(2)-(Ydim(2)-Ydim(1))/10)*ones(length(sigindices),1),'*','Color',[0,1,1],'FontSize',15,'HorizontalAlignment','center'); %display significance star
        end
    end
    
    % Display number of trials per average
    if showN
        Ydim = get(gca,'YLim');
        if showControl
            text(1,Ydim(1)+(Ydim(2)-Ydim(1))/10,sprintf('n=%d',rois(rindex).nTrials(s)),'HorizontalAlignment','center');
        end
        for s = 2:numStimuli
            text(s-~showControl,Ydim(1)+(Ydim(2)-Ydim(1))/10,sprintf('n=%d',rois(rindex).nTrials(s)),'HorizontalAlignment','center');
        end
    end
    
    % Display p-values
    if showPValues
        Ydim = get(gca,'YLim');
        for s = 1:numStimuli
            text(s,Ydim(1)+(Ydim(2)-Ydim(1))/20,sprintf('p=%.3f',rois(rindex).PValue(s)),'HorizontalAlignment','center');
        end
    end
    
    % Plot 0 line
    if showZero
        XLim = get(gca,'XLim');
        plot(XLim, [0 0], 'k--');
        temp = get(gca,'children');             % get ordering
        set(gca,'children',temp([2:end,1]));    % move line to back of plot
    end
    
    hold off
    
    % Save plot to PDF
    if saveOut && ismember(index, lastROI)
        drawnow;
        export_fig(hF, saveFile, '-append');
        close(hF);
    end
    
end %cycle ROIs


