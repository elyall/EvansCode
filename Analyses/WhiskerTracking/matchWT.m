function AnalysisInfo = matchWT(WTImageFiles, ExperimentFile, WTtimestamps)

saveOut = false;
saveFile = '';

directory = pwd;

%% Parse input arguments
if ~exist('WTImageFiles', 'var')
    [WTImageFiles, p] = uigetfile({'*.raw;*.tif;*.pgm'}, 'Select Whisker Tracking Image Files', directory, 'MultiSelect', 'on');
    if isnumeric(WTImageFiles)
        return
    elseif iscell(WTImageFiles)
            WTImageFiles = fullfile(p, WTImageFiles);
    elseif ischar(WTImageFiles)
        WTImageFiles = {fullfile(p, WTImageFiles)};
    end
elseif ischar(WTImageFiles)
    if isdir(WTImageFiles) % select all files in input directory
        p = WTImageFiles;
        WTImageFiles = dir(p);
        WTImageFiles = fullfile(p, {WTImageFiles(~[WTImageFiles(:).isdir]).name});
    else
        WTImageFiles = {WTImageFiles};
    end
end

if ~exist('ExperimentFile','var') || isempty(ExperimentFile)
    [ExperimentFile, p] = uigetfile({'*.mat;*.exp'},'Choose Experiment file',directory);
    if isnumeric(ExperimentFile)
        return
    end
    ExperimentFile = fullfile(p,ExperimentFile);
end

if saveOut && isempty(saveFile)
    saveFile = ExperimentFile;
end

%% Load in variables
load(ExperimentFile, 'Experiment', 'AnalysisInfo', 'InputNames', '-mat');
numTrials = size(AnalysisInfo, 1);

if ~exist('WTtimestamps', 'var')
    WTtimestamps = loadTimeStamps(WTImageFiles);
end

%% Load WT frame trigger
% [p,fn,~] = fileparts(ExperimentFile);
% DataInFile = fullfile(p,strcat(fn,'.bin')); %[strtok(StimFile, '.'),'_datain.bin'];
% if ~exist(DataInFile, 'file')
%     error('Cannot locate DataIn file for experiment: %s', DataInFile);
%     return
% end
% DataInFID = fopen(DataInFile, 'r');
% DataIn = fread(DataInFID, inf, Experiment.saving.dataPrecision);
% fclose(DataInFID);
% 
% % Select specific data
% nInputChannels = numel(InputNames);
% if nInputChannels > 1 % if more than one channel is saved
%     DataIn = reshape(DataIn, nInputChannels, numel(DataIn)/nInputChannels)';
% end
% DataIn = DataIn(:, find(strcmp(InputNames, 'I_WhiskerTracker')));
% 
% % Determine number of frames captured in each trial
% numFramesCaptured = nan(numTrials, 1);
% DataIn = DataIn - [0;DataIn(1:end-1)];
% for tindex = 1:numTrials
%     currentTriggers = DataIn(AnalysisInfo.ExpStimScans(1):AnalysisInfo.ExpStimScans(2));
%     numFramesCaptured(tindex) = sum(currentTriggers>0);
% end


%% Fix errant timer
n = 9;
numScans = size(WTtimestamps,1);
for sindex = 1:n:numScans
    lastscan = min(sindex+n-1,numScans);
    vals = unique(WTtimestamps(sindex:lastscan, 1));
    for index = 1:numel(vals)
        indices = WTtimestamps(sindex:lastscan, 1) == vals(index);
        if sum(indices) == 1
            WTtimestamps(sindex+find(indices)-1, 1) = WTtimestamps(sindex+find(indices)-2);
        end
    end
end

%% Determine breaks in whisker tracking frames
secondcounter = diff(WTtimestamps(:,1),1);
breaks = find(abs(secondcounter) > 1 & secondcounter > -129);
numFramesSaved = [breaks(1); diff(breaks)];
while numel(numFramesSaved) > numTrials
    [~,index] = min(numFramesSaved);
    if index~=1
        numFramesSaved(index-1) = numFramesSaved(index-1) + numFramesSaved(index);
    else
        numFramesSaved(index+1) = numFramesSaved(index+1) + numFramesSaved(index);
    end
    numFramesSaved(index) = [];
end

%% Assign frames to files (assumes first trial in WT and .exp are the same)

% Initialize columns in table
if ~any(strcmp(AnalysisInfo.Properties.VariableNames, 'WTFiles'))
    AnalysisInfo = [AnalysisInfo, table(cell(numTrials, 1), 'VariableNames', {'WTFiles'})];
end

% Parse files to different trials
framecounter = 1;
for tindex = 1:numel(numFramesSaved)
    AnalysisInfo.WTFiles{tindex} = WTImageFiles(framecounter:framecounter + numFramesSaved(tindex) - 1);
    framecounter = framecounter + numFramesSaved(tindex);
end


%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'AnalysisInfo', '-mat', '-v7.3');
    else
        save(saveFile, 'AnalysisInfo', '-mat', '-append');
    end
end