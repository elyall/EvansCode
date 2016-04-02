function [TrialIndex, RunIndex] = determineRunning(AnalysisInfo, frames, varargin)

TrialIndex = [1 inf];

MinThresh = -inf;
MeanThresh = 30;
MeanStdDevThresh = 1.3;
StdDevStdDevThresh = .8;

% Plots
verbose = false;
t = 1;

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'TrialIndex'
                TrialIndex = varargin{index+1};
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
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if (~exist('AnalysisInfo', 'var') && ~exist('frames', 'var')) || (isempty(AnalysisInfo) && isempty(frames))
    [ExperimentFile, p] = uigetfile('Select Experiment file:', directory);
    if ~ExperimentFile
        return
    end
    ExperimentFile = fullfile(p, ExperimentFile);
elseif ischar(AnalysisInfo)
    ExperimentFile = AnalysisInfo;
    clear AnalysisInfo
end


%% Load in variables
if ~exist('AnalysisInfo', 'var') || isempty(AnalysisInfo)
    load(ExperimentFile, 'AnalysisInfo', '-mat');
end
if ~exist('frames', 'var') || isempty(frames)
    load(ExperimentFile, 'frames', '-mat');
end

% Determine trials to analyze
if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(end-1)+1:size(AnalysisInfo, 1)];
end
numTrials = numel(TrialIndex);


%% Pull out running data
Data = cell(numTrials,1);
RunSpeed = zeros(numTrials,2);
for tindex = 1:numTrials
    if exist('frames', 'var') && ~isempty(frames)
        Data{tindex} = frames.RunningSpeed(AnalysisInfo.ExpStimFrames(TrialIndex(tindex),1):AnalysisInfo.ExpStimFrames(TrialIndex(tindex),2));
    else
        Data{tindex} = AnalysisInfo.meanRunningSpeed{TrialIndex(tindex)}(AnalysisInfo.TrialStimFrames(TrialIndex(tindex),1):AnalysisInfo.TrialStimFrames(TrialIndex(tindex),2));
    end
    RunSpeed(tindex,1) = mean(Data{tindex});
    RunSpeed(tindex,2) = std(Data{tindex});
end

% Convert Data to a matrix
D = cellfun(@numel,Data);
RunData = nan(max(D), numTrials);
for tindex = 1:numTrials
    RunData(1:D(tindex),tindex) = Data{tindex};
end


%% Determine trials in which mouse was running

% Initialize output
RunIndex = true(numTrials, 1);

% Plot initial data
if verbose
    hF = figure('NumberTitle','off','Name','Determining Running Speed');
    XLim = [min(RunSpeed)',max(RunSpeed)'];
    YLim = nan(3,2);
    first = true;
    plotdata('All Trials');
    first = false;
end

% Baseline min threshold
RunIndex(any(RunData < MinThresh)) = false;
if verbose; plotdata('After Thresholding Min'); end

% Baseline mean threshold
RunIndex(RunSpeed(:,1) < MeanThresh) = false;
if verbose; plotdata('After Thresholding Mean'); end

% Mean threshold
% while true
%     numTrials = nnz(RunIndex);
%     Indices = find(RunIndex);
%     Val = nan(numTrials,1);
%     for tindex = 1:numTrials
%         mu = mean(RunSpeed(Indices(setdiff(1:numTrials,tindex)),1));        % determine mean without current trial
%         sigma = std(RunSpeed(Indices(setdiff(1:numTrials,tindex)),1));      % determine std without current trial
%         Val(tindex) = (RunSpeed(Indices(tindex),1)-mu)/sigma;               % calculate zscore of current trial relative to other data
%     end
%     if any(Val < -MeanStdDevThresh)                                          % at least one outlier exists
%         [~,furthestIndex] = min(Val);                                       % determine largest outlier
%         RunIndex(Indices(furthestIndex)) = false;                           % remove largest outlier
%     else
%         break
%     end
%     if verbose
%         plotdata('Z-scoring Mean');
%     end
% end
RunIndex(RunSpeed(:,1) < mean(RunSpeed(RunIndex,1)) - MeanStdDevThresh*std(RunSpeed(RunIndex,1))) = false; % scott method
if verbose; plotdata('After Thresholding Mean Variance'); end

% Std dev threshold
% while true
%     numTrials = nnz(RunIndex);
%     Indices = find(RunIndex);
%     Val = nan(numTrials,1);
%     for tindex = 1:numTrials
%         mu = mean(RunSpeed(Indices(setdiff(1:numTrials,tindex)),2));        % determine mean without current trial
%         sigma = std(RunSpeed(Indices(setdiff(1:numTrials,tindex)),2));      % determine std without current trial
%         Val(tindex) = (RunSpeed(Indices(tindex),2)-mu)/sigma;               % calculate zscore of current trial relative to other data
%     end
%     if any(Val > StdDevStdDevThresh)                                        % at least one outlier exists
%         [~,furthestIndex] = max(Val);                                       % determine largest outlier
%         RunIndex(Indices(furthestIndex)) = false;                           % remove largest outlier
%     else
%         break
%     end
%     if verbose
%         plotdata('Z-scoring Std Dev');
%     end
% end
RunIndex(RunSpeed(:,2) > mean(RunSpeed(RunIndex,2)) + StdDevStdDevThresh*std(RunSpeed(RunIndex,2))) = false; % scott method
if verbose; plotdata('After Thresholding Std Dev Variance'); end

% Remove non-running trials from index
TrialIndex(~RunIndex) = [];


    function plotdata(Title)
        figure(hF);
        set(gcf,'Name',Title);
        
        subplot(1,3,1);
        plot(RunData(:,RunIndex));
        ylabel('Velocity (deg/s)');
        xlabel('Frame');
        title('Trial Speeds');
        axis tight;
        if first
            YLim(1,:) = get(gca,'YLim');
        else
            ylim(YLim(1,:));
        end
        
        subplot(1,3,2);
        hist(RunSpeed(RunIndex,1));
        ylabel('# of Trials');
        xlabel('Velocity (deg/s)');
        title('Trial Means');
        axis tight;
        xlim(XLim(1,:));
        if first
            YLim(2,:) = get(gca,'YLim');
        else
            ylim(YLim(2,:));
        end
        
        subplot(1,3,3);
        hist(RunSpeed(RunIndex,2));
        ylabel('# of Trials');
        xlabel('Std Dev (deg/s)');
        title('Trial Std Dev''s');
        axis tight;
        xlim(XLim(1,:));
        if first
            YLim(3,:) = get(gca,'YLim');
        else
            ylim(YLim(3,:));
        end
        
        drawnow;
        pause(t);
    end %plotdata

end %determineRunning
