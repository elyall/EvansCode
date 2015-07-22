function hA = plotTuningCurve2(rois, ROIindex, varargin)

saveToPDF = false;
filename = 'test.pdf'; %filename to save plots to

% Items to display
showFit = false;
showDataPoints = false;
showStimStars = false;
showN = false;
showPValues = false;

% Plot colors & display options
curveColor = 'r';
fitColor = 'g';
legendLocation = 'NorthWest'; % 'NorthWest' or 'NorthEast'
fontSize = 12;

% Constant variables
hA = [];
Title = [];

%% Check input arguments
if ~exist('rois', 'var') || isempty(rois)
    directory = cd;
    [rois, p] = uigetfile({'*.mat'},'Choose ROI file(s)', directory);
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
            legendLocation = varargin{index+1};
            index = index + 2;
        case 'fontSize'
            fontSize = varargin{index+1};
            index = index + 2;
        case 'axes'
            hA = varargin{index+1};
            index = index + 2;
        case 'Title'
            Title = varargin{index+1};
            index = index + 2;
        case 'saveToPDF'
            saveToPDF = true;
            index = index + 1;
        case 'filename'
            filename = varargin{index+1};
            index = index + 2;
        otherwise
            warning('Argument ''%s'' not recognized',varargin{index});
            index = index + 1;
    end
end

%% Load in data
if ischar(rois)
    ROIFile = rois;
    if saveToPDF && strcmp(filename, 'test.pdf')
        filename = strcat(ROIFile(1:end-3), 'pdf');
    end
    load(rois, 'ROIdata', '-mat');
    rois = ROIdata.rois;
end

%% Determine ROIs to plot
if ischar(ROIindex) && strcmp(ROIindex, 'all')
    ROIindex = 1:numel(rois);
end

%% Create axes
if isempty(hA)
    for rindex = 1:numel(ROIindex)
        figure();
        hA(rindex) = axes();
    end
end

%% Plot each specified ROI's tuning

numROIs = numel(ROIindex);
for index = 1:numROIs
        rindex = ROIindex(index);

        % Fix Label
        if isempty(rois(rindex).label)
            rois(rindex).label = {'none'};
        end
        
        % Select axes
        axes(hA{index}); %set axes as current axes
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
        errorbar(1,rois(rindex).curve(1),rois(rindex).StdError(1),curveColor,'LineStyle','-','LineWidth',2,'Marker','.','MarkerSize',20);  %plot control position
        errorbar(2:numStimuli,rois(rindex).curve(2:end),rois(rindex).StdError(2:end),curveColor,'LineStyle','-','LineWidth',2,'Marker','.','MarkerSize',20); %plot curve
        
        % Plot fit
        if showFit && numStimuli > 1
            yfit = feval(rois(rindex).Fit,1:0.001:numStimuli-1);
            plot(2:0.001:numStimuli, yfit + rois(rindex).offset, 'Color', fitColor, 'LineStyle', '-','LineWidth',1.25);
            Ydim = get(gca,'YLim');
            Xdim = get(gca,'XLim');
            switch legendLocation
                case 'NorthEast'
                    text(Xdim(2),Ydim(2),sprintf('r^2=%.2f\nFWHM=%.2f', rois(rindex).rsquare, rois(rindex).Coeff(3)),'Color',fitColor,'FontSize',fontSize,'HorizontalAlignment','right','VerticalAlignment','top');
                case 'NorthWest'
                    text(Xdim(1),Ydim(2),sprintf('r^2=%.2f\nFWHM=%.2f', rois(rindex).rsquare, rois(rindex).Coeff(3)),'Color',fitColor,'FontSize',fontSize,'HorizontalAlignment','left','VerticalAlignment','top');
            end
        end
        
        % Set axes labels
        set(gca,'XTick',1:numStimuli,'XTickLabel',[{'control'}; cellstr(num2str((1:numStimuli-1)'))]);
        xlabel('Position');
        ylabel('Average Stimulus-Evoked dF/F');
        xlim([0,numStimuli+1]);
        
        % Set Title
        if isempty(Title)
            title(sprintf('ROI: %d, Label: %s',rindex,rois(rindex).label{1}));
        else
            title(Title{rindex});
        end
        % title(sprintf('ROI: %d, Label: %s, File: %s',index,ROIdata.rois(rindex).label{1},ROIFiles{findex}));
        
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
        
        hold off
      
        % Save plot to PDF
        if saveToPDF
            drawnow;
            export_fig(figHandles{pindex}, filename, '-append');
            close(figHandles{pindex});
        end

end %cycle ROIs

% for pindex = 1:numFigs
%     figHandles{pindex};
%     legend('Pre', 'Post');
% end

