function [AnalysisInfo,frames,InputNames] = alignWTtriggers(ExperimentFile, varargin)

saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
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

if ~exist('ExperimentFile', 'var') || isempty(ExperimentFile) % nothing input
    [ExperimentFile,p] = uigetfile({'*.exp;*.mat'},'Choose Experiment File to process', directory);
    if isnumeric(ExperimentFile) % no file selected
        return
    end
    ExperimentFile = fullfile(p,ExperimentFile);
end

if saveOut && isempty(saveFile)
    saveFile = ExperimentFile;
end


%% Load in experiment data

% Load experiment info
load(ExperimentFile, 'DAQChannels', 'Experiment', 'TrialInfo', 'numDelayScans', '-mat');
% OutputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'O_')));
InputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
nInputChannels = numel(InputNames);
numTrials = numel(TrialInfo.StimID);

% Load binary DAQ data
[p,fn,~] = fileparts(ExperimentFile);
DataInFile = fullfile(p, strcat(fn,'.bin')); %[strtok(StimFile, '.'),'_datain.bin'];
DataInFID = fopen(DataInFile, 'r');
DataIn = fread(DataInFID, inf, Experiment.saving.dataPrecision);
fclose(DataInFID);
DataIn = reshape(DataIn, nInputChannels, numel(DataIn)/nInputChannels)';

% Load in table
load(ExperimentFile,'AnalysisInfo','frames','-mat');


%% Determine alignment of Whisker Tracking Frames
if any(ismember(InputNames, 'I_WhiskerTracker'))
    
    % Determine number of frames registered
    data = DataIn(:,strcmp(InputNames, 'I_WhiskerTracker'));
    data = data - [0;data(1:end-1)];
    data = data>0;
    
    frames.WhiskerTracking = nan(length(frames.Stimulus),1);
    for frameNumber = AnalysisInfo.ExpFrames(1,1):find(~isnan(frames.StartScan),1,'last')-1
        frames.WhiskerTracking(frameNumber) = sum(data(frames.StartScan(frameNumber):frames.StartScan(frameNumber+1)-1));
    end
    
    % Assign frames to specific trigger
    data = find(data);
    try
        AnalysisInfo = [AnalysisInfo, table(cell(numTrials,1),'VariableNames',{'WTStartScan'})];
    end
    for tindex = 1:numTrials
        AnalysisInfo.WTStartScan{tindex} = data(data>=AnalysisInfo.ExpScans(tindex,1) & data<=AnalysisInfo.ExpScans(tindex,2));
    end
    
end


%% Save output
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'AnalysisInfo', 'frames', '-mat', '-v7.3');
    else
        save(saveFile, 'AnalysisInfo', 'frames', '-mat', '-append');
    end
    fprintf('\tAnalysisInfo saved to: %s\n', saveFile);
end