function [TrialIndex, RunIndex, FileIndex] = determineRunning_Neuron2017(RunData, varargin)

TrialIndex = [1 inf]; %input data or data to load
StimID = [];
outFormat = '';

% Evan Method
% MinThresh = -inf;
% MeanThresh = 100;
% MeanStdDevThresh = 2;
% StdDevStdDevThresh = 2;

% Scott Method
MinThresh = -inf;
MeanThresh = 30;
MeanStdDevThresh = 1.3;
StdDevStdDevThresh = .8;

% Plots
verbose = false;
t = 0;

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

%% Mean and Std (Evan)

% % Std dev threshold
% while true
%     numTrials = nnz(RunIndex);
%     Indices = find(RunIndex);
%     Val = nan(numTrials,1);
%     for tindex = 1:numTrials
%         mu = mean(Metrics(Indices(setdiff(1:numTrials,tindex)),2));        % determine mean without current trial
%         sigma = std(Metrics(Indices(setdiff(1:numTrials,tindex)),2));      % determine std without current trial
%         Val(tindex) = (Metrics(Indices(tindex),2)-mu)/sigma;               % calculate zscore of current trial relative to other data
%     end
%     if any(Val > StdDevStdDevThresh)                                        % at least one outlier exists
%         [~,furthestIndex] = max(Val);                                       % determine largest outlier
%         RunIndex(Indices(furthestIndex)) = false;                           % remove largest outlier
%     else
%         break
%     end
% %     if verbose
% %         plotdata('Z-scoring Std Dev');
% %     end
% end
% if verbose; plotdata('After Thresholding Std Dev Variance'); end
% 
% % Mean threshold
% while true
%     numTrials = nnz(RunIndex);
%     Indices = find(RunIndex);
%     Val = nan(numTrials,1);
%     for tindex = 1:numTrials
%         mu = mean(Metrics(Indices(setdiff(1:numTrials,tindex)),1));        % determine mean without current trial
%         sigma = std(Metrics(Indices(setdiff(1:numTrials,tindex)),1));      % determine std without current trial
%         Val(tindex) = abs(Metrics(Indices(tindex),1)-mu)/sigma;            % calculate zscore of current trial relative to other data
%     end
%     if any(Val > MeanStdDevThresh)                                         % at least one outlier exists
%         [~,furthestIndex] = max(Val);                                      % determine largest outlier
%         RunIndex(Indices(furthestIndex)) = false;                          % remove largest outlier
%     else
%         break
%     end
% %     if verbose
% %         plotdata('Z-scoring Mean');
% %     end
% end
% if verbose; plotdata('After Thresholding Mean Variance'); end


%% Mean & Std (Scott)
RunIndex(Metrics(:,1) < mean(Metrics(RunIndex,1)) - MeanStdDevThresh*std(Metrics(RunIndex,1))) = false; % scott method
if verbose; plotdata('After Thresholding Mean Variance'); end
RunIndex(Metrics(:,2) > mean(Metrics(RunIndex,2)) + StdDevStdDevThresh*std(Metrics(RunIndex,2))) = false; % scott method
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
        temp(cellfun(@isempty,temp)) = {0};
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
        pause(t);
    end %plotdata

end %determineRunning
