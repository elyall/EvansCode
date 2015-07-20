function WTdata = processWT(WTImageFiles, ExperimentFile)

saveOut = false;
saveFile = '';


%% Parse input arguments
if ~exist('WTImageFiles', 'var')
    [WTImageFiles, p] = uigetfile({'*.raw;*.tif;*.pgm'}, 'Select Whisker Tracking Image Files', directory, 'MultiSelect', 'on');
    if isnumeric(WTImageFiles)
        return
    elseif iscell(WTImageFiles)
        for findex = 1:numel(WTImageFiles)
            WTImageFiles{findex} = fullfile(p, WTImageFiles{findex});
        end
    elseif ischar(WTImageFiles)
        WTImageFiles = {fullfile(p, WTImageFiles)};
    end
elseif ischar(WTImageFiles)
    WTImageFiles = {WTImageFiles};
end

if ~exist('ExperimentFile','var') || isempty(ExperimentFile)
    directory = CanalSettings('ExperimentDirectory');
    [ExperimentFile, p] = uigetfile({'*.mat;*.exp'},'Choose Experiment file',directory);
    if isnumeric(ExperimentFile)
        return
    end
    ExperimentFile = fullfile(p,ExperimentFile);
end

%% Load in variables
load(ExperimentFile, 'Experiment', 'AnalysisInfo', 'InputNames');

%% Load WT frame trigger
[p,fn,~] = fileparts(ExperimentFile);
DataInFile = fullfile(p,strcat(fn,'.bin')); %[strtok(StimFile, '.'),'_datain.bin'];
if ~exist(DataInFile, 'file')
    error('Cannot locate DataIn file for experiment: %s', DataInFile);
    return
end
DataInFID = fopen(DataInFile, 'r');
DataIn = fread(DataInFID, inf, Experiment.saving.dataPrecision);
fclose(DataInFID);

% Select specific data
nInputChannels = numel(InputNames);
if nInputChannels > 1 % if more than one channel is saved
    DataIn = reshape(DataIn, nInputChannels, numel(DataIn)/nInputChannels)';
end
DataIn = DataIn(:,find(strcmp(InputNames, 'I_WhiskerTracker')));

%% Determine which files occurred during each trial
DataIn = DataIn - [0;DataIn(1:end-1)];
numFrames = sum(DataIn>0);
numFrames = numel(WTImageFiles)
% for tindex = 1:size(AnalysisInfo, 1)
%     currentTriggers = DataIn(AnalysisInfo.ExpStimScans(1):AnalysisInfo.ExpStimScans(2));
%     currentTriggers = currentTriggers - [0;currentTriggers(1:end-1)];
%     numFrames = sum(currentTriggers>0);
%     
%     
% end
