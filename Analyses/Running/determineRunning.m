function [TrialIndex, RunIndex, FileIndex] = determineRunning(RunData, varargin)


% Default parameters that can be adjusted
MinThresh = -inf;           % scalar that results in any trial with an instantaneous velocity below this value being thrown out (set to -inf to not throw anything out)
MeanThresh = 100;           % scalar that results in any trial with a mean velocity below this value will be thrown out (set to -inf to not throw anything out)
MeanMethod = 'medianRule';  % 'tukey','evan', or 'medianRule', determines the method used to find outliers within the remaining distribution of trial mean velocities
StdMethod = 'medianRule';   % 'tukey','evan', or 'medianRule', determines the method used to find outliers within the remaining distribution of trial std velocities
TrialIndex = [1 inf];       % vector specifying the indices of the input data or data to load
outFormat = '';             % 'numeric' or 'cell' specifying whether the output is a vector of trial indices or a cell containing said vector

% Plotting
verbose = false;            % display plots
StimID = [];                % vector specifying the stim indices of each trial (only used for plotting)
PauseTime = 2;              % time to pause after displaying each plot

% Evan's method parametersw
numSTDsOutlier = 4;     % inf if no threshold
minNumTrials = 5;       % -inf if no minimum (will error when gets to 2 trials)

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'TrialIndex'
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'StimID'
                StimID = varargin{index+1};
                index = index + 2;
            case 'outFormat'
                outFormat = varargin{index+1};
                index = index + 2;
            case 'MinThresh'
                MinThresh = varargin{index+1};
                index = index + 2;
            case 'MeanThresh'
                MeanThresh = varargin{index+1};
                index = index + 2;
            case 'MeanStdDevThresh'
                MeanStdDevThresh = varargin{index+1};
                index = index + 2;
            case 'StdDevStdDevThresh'
                StdDevStdDevThresh = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
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

if ~exist('RunData','var') || isempty(RunData)
    [RunData, p] = uigetfile({'*.exp'},'Choose Experiment file',directory,'MultiSelect','on');
    if isnumeric(RunData)
        return
    end
    RunData = fullfile(p,RunData);
end
if ~iscell(RunData)
    RunData = {RunData};
    if isempty(outFormat)
        outFormat = 'numeric';
    end
elseif isempty(outFormat)
    outFormat = 'cell';
end
numFiles = numel(RunData);

if ~iscell(TrialIndex)
    TrialIndex = {TrialIndex};
end
if numel(TrialIndex) == 1 && numFiles > 1
    TrialIndex = repmat(TrialIndex,numFiles,1);
end


%% Load in variables
if iscellstr(RunData)
    StimID = [];
    for findex = 1:numFiles
        [RunData{findex},TrialIndex{findex},temp] = gatherRunData(RunData{findex},[],'type','mat','TrialIndex',TrialIndex{findex});
        if verbose
            if findex ~= 1
                temp = temp + max(StimID) + 1;
            end
            StimID = [StimID;temp];
        end
    end
end

nT = cellfun(@numel,TrialIndex);
numTrials = sum(nT);

% Format data
FileIndex = cell(numFiles,1);
sz = cellfun(@size,RunData,'UniformOutput',false);
sz = reshape([sz{:}],2,numFiles)';
m = max(sz(:,2));
for findex = 1:numFiles
    RunData{findex} = [RunData{findex},nan(sz(findex,1),m-sz(findex,2))]';
    FileIndex{findex} = findex*ones(nT(findex),1)';
end
RunData = [RunData{:}]'; % merge run data to one matrix
FileIndex = [FileIndex{:}]';

if verbose && isempty(StimID)
    StimID = ones(numTrials,1);
end


%% Pull out running data
Metrics = [nanmean(RunData,2), nanstd(RunData,[],2)];


%% Determine trials in which mouse was running

% Initialize output
RunIndex = true(numTrials, 1);

% Plot initial data
if verbose
    hF = figure('NumberTitle','off','Name','Determining Running Speed','Position',[50, 50, 1400, 800]);
    StimIDs = unique(StimID);
    numStims = numel(StimIDs);
    Colors = lines(numStims);
    XLim = nan(5,2);
    YLim = nan(5,2);
    first = true;
    plotdata('All Trials');
    first = false;
end

%% Baseline min threshold
RunIndex(any(RunData < MinThresh,2)) = false;
if verbose; plotdata('After Thresholding Min'); end

%% Baseline mean threshold
RunIndex(Metrics(:,1) < MeanThresh) = false;
if verbose; plotdata('After Thresholding Mean'); end

%% Determine mean outliers
outliers = determineOutliers(Metrics(RunIndex,1), 'type', MeanMethod, 'numSTDsOutlier', numSTDsOutlier, 'minNumTrials', minNumTrials);
temp = find(RunIndex);
RunIndex(temp(outliers)) = false;
if verbose; plotdata('After Thresholding Mean Variance'); end

%% Determine std outliers
outliers = determineOutliers(Metrics(RunIndex,2), 'type', StdMethod, 'numSTDsOutlier', numSTDsOutlier, 'minNumTrials', minNumTrials);
temp = find(RunIndex);
RunIndex(temp(outliers)) = false;
if verbose; plotdata('After Thresholding Std Dev Variance'); end

%% Remove non-running trials from index
for findex = 1:numFiles
    TrialIndex{findex}(~RunIndex(FileIndex==findex)) = []; % remove non-running trials from input TrialIndex
end
if strcmp(outFormat,'numeric')
    TrialIndex = [TrialIndex{:}];
end


%% Plot Data
    function plotdata(Title)
        figure(hF);
        set(gcf,'Name',Title);
        
        % Raw velocities for each trial
        subplot(3,3,[1,4,7]);
        plot(RunData(~RunIndex,:)','Color',[.9,.9,.9]); hold on;
        for sindex = 1:numStims
            goodIndex = all([RunIndex,StimID==StimIDs(sindex)],2);
            plot(RunData(goodIndex,:)','Color',Colors(sindex,:));
        end
        hold off;
        ylabel('Velocity (deg/s)');
        xlabel('Frame');
        title(sprintf('Trial Speeds (%d/%d:%.f%%)',nnz(RunIndex),numTrials,nnz(RunIndex)/numTrials*100));
        axis tight;
        if first
            XLim(1,:) = get(gca,'XLim');
            YLim(1,:) = get(gca,'YLim');
        else
            ylim(YLim(1,:));
            xlim(XLim(1,:));
        end
        
        % Trial means
        subplot(3,3,2);
        histogram(Metrics(~RunIndex,1),'BinWidth',50); hold on;
        histogram(Metrics(RunIndex,1),'BinWidth',50); hold off;
        ylabel('# of Trials');
        xlabel('Mean Velocity (deg/s)');
        title('Trial Mean Velocities');
        axis tight;
        if first
            XLim(2,:) = get(gca,'XLim');
            YLim(2,:) = get(gca,'YLim');
        else
            xlim(XLim(2,:));
            ylim(YLim(2,:));
        end    
        
        % Trial std devs
        subplot(3,3,3);
        histogram(Metrics(~RunIndex,2),'BinWidth',20); hold on;
        histogram(Metrics(RunIndex,2),'BinWidth',20); hold off;
        ylabel('# of Trials');
        xlabel('Std Dev (deg/s)');
        title('Trial Standard Deviations');
        axis tight;
        if first
            XLim(3,:) = get(gca,'XLim');
            YLim(3,:) = get(gca,'YLim');
        else
            xlim(XLim(3,:));
            ylim(YLim(3,:));
        end
        
        % Trial means over time
        subplot(3,3,[5,6]);
        plot(find(~RunIndex),Metrics(~RunIndex,1),'b.'); hold on;
        plot(find(RunIndex),Metrics(RunIndex,1),'r.','MarkerSize',5);
        plot(find(RunIndex),smooth(find(RunIndex),Metrics(RunIndex,1)),'g--'); 
        coeffs = polyfit(find(RunIndex), Metrics(RunIndex,1), 1);
        fittedY = polyval(coeffs, find(RunIndex));
        plot(find(RunIndex),fittedY,'k');
        hold off;
        ylabel('Mean Velocity (deg/s)');
        xlabel('Trial #');
        title('Mean Velocity vs. Time');
        axis tight;
        if first
            XLim(4,:) = get(gca,'XLim');
            YLim(4,:) = get(gca,'YLim');
        else
            xlim(XLim(4,:));
            ylim(YLim(4,:));
        end
        
        % Trial means per stim
        subplot(3,3,[8,9]);
        temp = cell(numStims,1);
        for sindex = 1:numStims
            temp{sindex} = Metrics(all([RunIndex,StimID==StimIDs(sindex)],2),1);
        end
        temp(cellfun(@isempty,temp)) = {0}; % fix empty cells
        violin(temp'); grid on; hold on;
        for sindex = 1:numStims
            badIndex = all([~RunIndex,StimID==StimIDs(sindex)],2);
            plot((sindex-.1)*ones(nnz(badIndex),1),Metrics(badIndex,1),'b.');
            goodIndex = all([RunIndex,StimID==StimIDs(sindex)],2);
            plot((sindex+.1)*ones(nnz(goodIndex),1),Metrics(goodIndex,1),'r.');
            if ~first
                text(sindex,YLim(5,1)+.1*range(YLim(5,:)),sprintf('%d/%d',nnz(goodIndex),(nnz(goodIndex)+nnz(badIndex))),'HorizontalAlignment','Center');
            end
        end
        hold off;
        ylabel('Mean Velocity (deg/s)');
        xlabel('Stimulus');
        title('Mean Velocity vs. Stimulus');
        axis tight;
        if first
            XLim(5,:) = get(gca,'XLim') + [-.2,.2];
            YLim(5,:) = get(gca,'YLim');
        else
            xlim(XLim(5,:));
            ylim(YLim(5,:)*.8);
        end     
        
        drawnow;
        pause(PauseTime);
    end %plotdata

end %determineRunning
