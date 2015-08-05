function [hA, hF] = plotTuningCurve(rois, ROIindex, varargin)

saveOut = false;
saveFile = ''; % filename to save plots to

% Items to display
showFit = false;
showDataPoints = false;
showStimStars = false;
showN = false;
showPValues = false;

% Plot colors & display options
curveColor = 'r';
fitColor = 'g';
fitLegendLocation = 'NorthWest'; % 'NorthWest' or 'NorthEast'
fontSize = 12;

% Constant variables
FigureIndex = [];
hF = [];
AxesIndex = [];
hA = [];
SubplotDim = [];
Title = [];
YLim = [];
Legend = {'Pre Control', 'Pre', 'Post Control', 'Post'};
LegendLocation = 'NorthEastOutside';

directory = cd;

%% Parse input arguments
if ~exist('rois', 'var') || isempty(rois)
    directory = cd;
    [rois, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory);
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
            case 'FigureIndex'
                FigureIndex = varargin{index+1};
                index = index + 2;
            case 'figures'
                hF = varargin{index+1};
                index = index + 2;
            case 'AxesIndex'
                AxesIndex = varargin{index+1};
                index = index + 2;
            case 'axes'
                hA = varargin{index+1};
                index = index + 2;
            case 'SubplotDim'
                SubplotDim = varargin{index+1};
                index = index + 2;
            case 'Title'
                Title = varargin{index+1};
                index = index + 2;
            case 'YLim'
                YLim = varargin{index+1};
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


%% Determine number of figures
if isempty(FigureIndex)
    FigureIndex = 1:numROIs;
    hF = nan(numROIs, 1);
elseif isempty(hF)
    hF = nan(max(FigureIndex), 1);
end
numFigs = numel(hF);


%% Determine number of axes
if isempty(AxesIndex)
    AxesIndex = ones(numROIs, 1);
    hA = cell(numFigs, 1);
    [hA{:}] = deal(nan);
end

% Determine subplot configuration for each figure
numAxesPerFig = zeros(numFigs, 1);
if isempty(SubplotDim)
    SubplotDim = ones(numFigs, 2);
    popSubplot = true;
else
    popSubplot = false;
end
for findex = 1:numFigs
    numAxesPerFig(findex) = max(AxesIndex(FigureIndex==findex));
    if popSubplot
       SubplotDim(findex, 2) = numAxesPerFig(findex);
    end
end

% Initialize handles container
if isempty(hA)
    hA = cell(numFigs, 1);
    for findex = 1:numFigs
        hA{findex} = nan(numAxesPerFig(findex));
    end
end

% Initialize titles container
if isempty(Title)
    Title = cell(numFigs, 1);
    for findex = 1:numFigs
        Title{findex} = cell(1,numAxesPerFig(findex));
    end
end

% Determine which ROI is first & last for each figure
[~, lastROI] = unique(FigureIndex, 'last');


%% Determine plotting colors
if ischar(curveColor)
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
        if isempty(rois(rindex).label)
            rois(rindex).label = {'none'};
        end
        
        % Create figure
        if isnan(hF(FigureIndex(index)))
            hF(FigureIndex(index)) = figure();
        else
            figure(hF(FigureIndex(index)));
        end
        
        % Select axes
        hA{FigureIndex(index)}(AxesIndex(index)) = subplot(SubplotDim(FigureIndex(index), 1), SubplotDim(FigureIndex(index), 2), AxesIndex(index)); %set axes as current axes
        hold on
        
        % Plot each trial's average dF/F for each stimulus
        if showDataPoints
            plot(ones(rois(rindex).nTrials(1),1), rois(rindex).Raw{1}, 'k.') %plot raw data points for control position
            for s = 2:rois(rindex).nstim
                plot((s)*ones(rois(rindex).nTrials(s),1), rois(rindex).Raw{s}, 'k.') %plot raw data points for all stimuli
            end
        end
        
        % Plot tuning curves
        numStimuli = numel(rois(rindex).curve);
        errorbar(1,rois(rindex).curve(1),rois(rindex).StdError(1),curveColor{index},'LineStyle','-','LineWidth',2,'Marker','.','MarkerSize',20);  %plot control position
        errorbar(2:numStimuli,rois(rindex).curve(2:end),rois(rindex).StdError(2:end),curveColor{index},'LineStyle','-','LineWidth',2,'Marker','.','MarkerSize',20); %plot curve
        
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
        set(gca,'XTick',1:numStimuli,'XTickLabel',[{'control'}; cellstr(num2str((1:numStimuli-1)'))]);
        xlabel('Position');
        ylabel('Average Stimulus-Evoked dF/F');
        xlim([0,numStimuli+1]);
        if ~isempty(YLim)
            ylim(YLim);
        end
        
        % Set Title
        if isempty(Title{FigureIndex(index)}{AxesIndex(index)})
            Title{FigureIndex(index)}{AxesIndex(index)} = sprintf('ROI: %d, Label: %s',rindex,rois(rindex).label{1});
        end
        title(Title{FigureIndex(index)}{AxesIndex(index)});
        
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
            for s = 1:numStimuli
                text(s,Ydim(1)+(Ydim(2)-Ydim(1))/10,sprintf('n=%d',rois(rindex).nTrials(s)),'HorizontalAlignment','center');
            end
        end
        
        % Display p-values
        if showPValues
            Ydim = get(gca,'YLim');
            for s = 1:numStimuli
                text(s,Ydim(1)+(Ydim(2)-Ydim(1))/20,sprintf('p=%.3f',rois(rindex).PValue(s)),'HorizontalAlignment','center');
            end
        end
        
        % Place legend
        if ~isempty(Legend) && ismember(index, lastROI)
            legend(Legend, 'Location', LegendLocation);
        end
        
        hold off
        
        % Save plot to PDF
        if saveOut && ismember(index, lastROI)
            drawnow;
            export_fig(hF(FigureIndex(index)), saveFile, '-append');
            close(hF(FigureIndex(index)));
        end

end %cycle ROIs


