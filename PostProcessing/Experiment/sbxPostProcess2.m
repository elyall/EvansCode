function [AnalysisInfo,frames,InputNames,Config] = sbxPostProcess2(StimFile, ImageFile, saveout)
% Use to be 'FormatExperiment': Format data to be analyzed

saveInputData = false;
numBaselineFrames = 'all'; % # of frames to record prior to stimulus, or 'all' to record every earlier frame in trial

%% Initialization
narginchk(0, 3);
if ~exist('StimFile', 'var') || isempty(StimFile)
    directory = CanalSettings('ExperimentDirectory');
    [StimFile,p] = uigetfile({'*.mat'}, 'Choose Stimulation file', directory);
    if isnumeric(StimFile)
        return
    end
    directory = p;
    StimFile = fullfile(p,StimFile);
end

if ~exist('ImageFile', 'var') || isempty(ImageFile)
    directory = CanalSettings('DataDirectory');
    [ImageFile,p] = uigetfile({'*.sbx'}, 'Choose sbx file', directory);
    if isnumeric(ImageFile)
        return
    end
    ImageFile = fullfile(p,ImageFile);
end

if nargin<3
    saveout = false;
    SaveFile = StimFile;
elseif ischar(saveout)
    SaveFile = savout;
    saveout = true;
elseif islogical(saveout) && saveout == true
    SaveFile = StimFile;
end

%% Determine DataIn File
[p,fn,~] = fileparts(StimFile);
DataInFile = fullfile(p,strcat(fn,'.bin')); %[strtok(StimFile, '.'),'_datain.bin'];

%% Determine image info
[~,~,ext] = fileparts(ImageFile);
switch ext
    case {'.sbx','.imgs'}
        InfoFile = sbxIdentifyFiles(ImageFile);
        InfoFile = InfoFile{1};
    case '.mat'
        InfoFile = ImageFile;
        ImageFile = sbxIdentifyFiles(InfoFile);
        ImageFile = ImageFile{1};
end
Config = parseSbxHeader(ImageFile);

%% Load Data

% Load Data Files
load(StimFile, 'DAQChannels', 'Experiment', 'TrialInfo', '-mat');
info = fixTimeStamps(InfoFile); %for data collected w/ Scanbox v1.2
% load(InfoFile, 'info');
if ischar(DataInFile)
    DataInFID = fopen(DataInFile, 'r');
    DataIn = fread(DataInFID, inf, Experiment.saving.dataPrecision);
    fclose(DataInFID);
end


% Format DataIn
InputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
nInputChannels = numel(InputNames);
if ischar(DataInFile) && nInputChannels > 1 % if more than one channel is saved
    DataIn = reshape(DataIn, nInputChannels, numel(DataIn)/nInputChannels)';
end


%% Creat Table With Row For Each Frame
scansPerFrame = Experiment.params.samplingFrequency/Config.FrameRate;
if saveout
    VariablesToSave = {'AnalysisInfo', 'scansPerFrame', 'InputNames', 'frames'};
end

nTrials = numel(TrialInfo.StimID);
nTrials_imaging = numel(info.frame(info.event_id==1))/4;
nTrialsMissed = nTrials - nTrials_imaging;
IndexStim = 1:nTrials;
IndexImaging = zeros(nTrials,1);
IndexImaging(IndexStim(nTrialsMissed+1:end)) = 1:nTrials_imaging;

AnalysisInfo = table(zeros(nTrials,1),zeros(nTrials,1),zeros(nTrials,1),cell(nTrials,1),zeros(nTrials,1),zeros(nTrials,1),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),...
    'VariableNames', {'StimID', 'TrialIndex', 'ImgIndex', 'ImgFilename', 'nFrames', 'nScans', 'ExpStimFrames', 'ExpStimScans', 'StimFrameLines', 'TrialStimFrames', 'TrialStimScans', 'ExpFrames', 'ExpScans'});

frames = struct('Stimulus', nan(Config.Frames,1), 'Trial', nan(Config.Frames,1));

if any(ismember(InputNames, 'I_RunWheelA'))
    AnalysisInfo = [AnalysisInfo, table(cell(nTrials,1),'VariableNames',{'meanRunningSpeed'})];
    frames.RunningSpeed = nan(Config.Frames,1);
end
if saveInputData   
    AnalysisInfo = [AnalysisInfo, table(cell(nTrials,1),'VariableNames',{'InputData'})];
    VariablesToSave = [VariablesToSave, {'InputNames'}];
end
if ischar(numBaselineFrames) || numBaselineFrames > 0
    AnalysisInfo = [AnalysisInfo, table(zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),'VariableNames',{'ExpBaselineFrames','ExpBaselineScans','TrialBaselineFrames','TrialBaselineScans'})];
end

scanboxFrames = info.frame(info.event_id==1);
scanboxLines = info.line(info.event_id==1);
for t = 1:nTrials
    AnalysisInfo.StimID(t) = TrialInfo.StimID(t); % ID # of stimulus presented
    AnalysisInfo.TrialIndex(t) = t; % index of trial presented relative to 'TrialInfo'
    AnalysisInfo.ImgIndex(t) = IndexImaging(t); % index of trial presented relative to what the images captured
    AnalysisInfo.ImgFilename{t} = Config.FullFilename(1); % filename of imaging file
    AnalysisInfo.nScans(t) = numel(Experiment.Stimulus); % number of scans in trial
    if t ~= 1 % determine index of scan that occurred last in previous trial
        lastscan = sum(AnalysisInfo.nScans(1:t-1));
    else % first trial
        lastscan = 0;
    end
    AnalysisInfo.ExpScans(t,:) = [lastscan + 1, lastscan + AnalysisInfo.nScans(t)]; % indices of scans for current trial relative to entire experiment
    if t > nTrialsMissed
        AnalysisInfo.ExpStimFrames(t,:) = scanboxFrames(4*(IndexImaging(t)-1)+2:4*IndexImaging(t)-1)'; % indices of first and last frame in which the stimulus occurred relative to entire experiment
        AnalysisInfo.StimFrameLines(t,:) = scanboxLines(4*(IndexImaging(t)-1)+2:4*IndexImaging(t)-1)'; % indices of first and last line in which the stimulus occurred relative to each frame
        %AnalysisInfo.ExpStimFrames(t,:) = info.frame(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % indices of first and last frame in which the stimulus occurred relative to entire experiment
        %AnalysisInfo.StimFrameLines(t,:) = info.line(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % indices of first and last line in which the stimulus occurred relative to each frame
        %AnalysisInfo.EventIDs(t,:) = info.event_id(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % ID # of TTL port in ScanBox that first and last trigger came in on
        frames.Stimulus(AnalysisInfo.ExpStimFrames(t,1):AnalysisInfo.ExpStimFrames(t,2)) = AnalysisInfo.StimID(t);
        AnalysisInfo.TrialStimScans(t,:) = [find(Experiment.Stimulus,1),find(Experiment.Stimulus,1,'last')]; % indices of first and last scan in which the stimulus occurred relative to current trial
        AnalysisInfo.TrialStimFrames(t,:) = [floor(AnalysisInfo.TrialStimScans(t,1)/scansPerFrame),ceil(AnalysisInfo.TrialStimScans(t,2)/scansPerFrame)]; % indices of first and last frame in which the stimulus occurred relative to current trial
        AnalysisInfo.ExpFrames(t,:) = [AnalysisInfo.ExpStimFrames(t,1) - AnalysisInfo.TrialStimFrames(t,1) + 1, min(AnalysisInfo.ExpStimFrames(t,1) - AnalysisInfo.TrialStimFrames(t,1) + ceil(numel(Experiment.Stimulus)/scansPerFrame),Config.Frames)]; % indices of frames for current trial relative to entire experiment
        frames.Trial(AnalysisInfo.ExpFrames(t,1):AnalysisInfo.ExpFrames(t,2)) = t;
        AnalysisInfo.nFrames(t) = AnalysisInfo.ExpFrames(t,2) - AnalysisInfo.ExpFrames(t,1) + 1;
        AnalysisInfo.ExpStimScans(t,:) = lastscan + AnalysisInfo.TrialStimScans(t,:);
    end
    if saveInputData
        AnalysisInfo.InputData{t} = DataIn(AnalysisInfo.ExpScans(t,1):min(end,AnalysisInfo.ExpScans(t,2)),:);
    end
    if any(ismember(InputNames, 'I_RunWheelA'))
        AnalysisInfo.meanRunningSpeed{t} = zeros(AnalysisInfo.nFrames(t),1);
        data = DataIn(AnalysisInfo.ExpScans(t,1):min(end,AnalysisInfo.ExpScans(t,2)),find(strcmp(InputNames, 'I_RunWheelB')));
        dataB = DataIn(AnalysisInfo.ExpScans(t,1):min(end,AnalysisInfo.ExpScans(t,2)),find(strcmp(InputNames, 'I_RunWheelA')));
        data = (data - [0; data(1:end-1)]);
        data(data<0) = 0;
        dataB(~data) = 0;
        data(logical(dataB)) = -1;
        if ~isempty(data)
            data = calcRunningSpeed(data, Experiment.params.samplingFrequency, 10, date);
            dataPointsPerFrame = round(scansPerFrame/10);
            for frameNumber = 1:AnalysisInfo.nFrames(t)
                AnalysisInfo.meanRunningSpeed{t}(frameNumber) = mean(data(dataPointsPerFrame*(frameNumber-1)+1:min(end,dataPointsPerFrame*frameNumber)));
                frames.RunningSpeed(AnalysisInfo.ExpFrames(t,1)+frameNumber-1) = AnalysisInfo.meanRunningSpeed{t}(frameNumber);
            end
        else
            for frameNumber = 1:AnalysisInfo.nFrames(t)
                 AnalysisInfo.meanRunningSpeed{t}(frameNumber) = NaN;
                 frames.RunningSpeed(AnalysisInfo.ExpFrames(t,1)+frameNumber-1) = AnalysisInfo.meanRunningSpeed{t}(frameNumber);
            end
        end
    end
    if ischar(numBaselineFrames) || numBaselineFrames > 0
        if ischar(numBaselineFrames)
            AnalysisInfo.ExpBaselineFrames(t,:) = [AnalysisInfo.ExpFrames(t,1), AnalysisInfo.ExpStimFrames(t,1)-1];
            AnalysisInfo.ExpBaselineScans(t,:) = [AnalysisInfo.ExpScans(t,1), AnalysisInfo.ExpStimScans(t,1) - 1];
            AnalysisInfo.TrialBaselineFrames(t,:) = [1, AnalysisInfo.TrialStimFrames(t,1) - 1];
            AnalysisInfo.TrialBaselineScans(t,:) = [1, AnalysisInfo.TrialStimScans(t,1) - 1];
        else
            AnalysisInfo.ExpBaselineFrames(t,:) = [max(1,AnalysisInfo.ExpStimFrames(t,1) - numBaselineFrames), AnalysisInfo.ExpStimFrames(t,1)-1];
            AnalysisInfo.ExpBaselineScans(t,:) = [max(1,AnalysisInfo.ExpStimScans(t,1) - ceil(numBaselineFrames*scansPerFrame)), AnalysisInfo.ExpStimScans(t,1) - 1];
            AnalysisInfo.TrialBaselineFrames(t,:) = [max(1,AnalysisInfo.TrialStimFrames(t,1) - numBaselineFrames), AnalysisInfo.TrialStimFrames(t,1) - 1];
            AnalysisInfo.TrialBaselineScans(t,:) = [max(1,AnalysisInfo.TrialStimScans(t,1) - ceil(numBaselineFrames*scansPerFrame)), AnalysisInfo.TrialStimScans(t,1) - 1];
        end
    end
end

% frames = struct2table(frames);

if saveout
    save(SaveFile, VariablesToSave{:}, '-mat', '-append');
    fprintf('AnalysisInfo saved to: %s\n', SaveFile);
end

