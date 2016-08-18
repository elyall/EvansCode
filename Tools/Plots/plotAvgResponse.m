function hA = plotAvgResponse(ROIdata, ROIindex, varargin)

hA = [];

TrialIndex = [];
StimIDs = [];
FrameIndex = [];

Title = '';
LabelAxes = false;

frameRate = 15.45;
LineWidth = 1;
LineStyle = '-';
Colors = [];

showLegend = false;

%% Parse input arguments
if ~exist('ROIdata', 'var') || isempty(ROIdata)
    [ROIdata, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if ~exist('ROIindex','var') || isempty(ROIindex)
    ROIindex = 1;
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Axes','axes','hA'}
                hA = varargin{index+1};
                index = index + 2;
            case 'TrialIndex'
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'FrameIndex'
                FrameIndex = varargin{index+1};
                index = index + 2;
            case 'StimIDs'
                StimIDs = varargin{index+1};
                index = index + 2;
            case 'LabelAxes'
                LabelAxes = true;
                index = index + 1;
            case 'Title'
                Title = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
                index = index + 2;
            case 'LineWidth'
                LineWidth = varargin{index+1};
                index = index + 2;
            case 'LineStyle'
                LineStyle = varargin{index+1};
                index = index + 2;
            case 'Colors'
                Colors = varargin{index+1};
                index = index + 2;
            case 'showLegend'
                showLegend = true;
                index = index + 1;
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
if ischar(ROIdata)
    load(ROIdata, 'ROIdata', '-mat');
end

%% Determine specifics

% determine trials
if isempty(TrialIndex)
    TrialIndex = ROIdata.DataInfo.TrialIndex;
end
TrialIndex = ismember(ROIdata.DataInfo.TrialIndex', TrialIndex);

% determine stimuli
if isempty(StimIDs)
    StimIDs = unique(ROIdata.DataInfo.StimID(TrialIndex));
elseif isrow(StimIDs)
    StimIDs = StimIDs';
end
numStims = numel(StimIDs);

% determine colors
if isempty(Colors)
    Colors = lines(numStims);
end

% determine frames
if isempty(FrameIndex)
    FrameIndex = [1 inf];
end
if FrameIndex(1) <= 1
    FrameIndex(1) = 1;
    xshift = 0;
else
    xshift = FrameIndex(1) - 1;
end
if FrameIndex(2) > size(ROIdata.rois(ROIindex).dFoF,2)
    FrameIndex(2) = size(ROIdata.rois(ROIindex).dFoF,2);
end
    

%% Plot data

% Create axis
if isempty(hA)
    figure;
    hA = axes;
else
    axes(hA);
end
hold on;

% Plot data
for sindex = 1:numStims
    data = mean(ROIdata.rois(ROIindex).dFoF(TrialIndex&ROIdata.DataInfo.StimID==StimIDs(sindex),FrameIndex(1):FrameIndex(2)));
    plot(0:1/frameRate:range(FrameIndex)/frameRate,data,'Color',Colors(sindex,:),'LineStyle',LineStyle,'LineWidth',LineWidth);
end
axis tight;
xlim([0,range(FrameIndex)/frameRate]);

% Plot stim lines
YLim = get(gca,'ylim');
f = ((ROIdata.DataInfo.numFramesBefore + .5 - xshift)-1)/frameRate;
plot([f,f],YLim,'k--');
l = ((ROIdata.DataInfo.numFramesBefore + mode(ROIdata.DataInfo.numStimFrames) + .5 - xshift)-1)/frameRate;
plot([l,l],YLim,'k--');
hold off;

% Label axes
if LabelAxes
    xlabel('Time (s)');
    ylabel('dF/F');
end

% Display title
if ~isempty(Title)
    title(Title);
end

% Display legend
if showLegend
    legend(cellstr(num2str(StimIDs)),'Location','NorthEastOutside');
end


