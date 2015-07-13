function [AnalysisInfo,frames,InputNames,ImageInfo] = sbxPostProcess(StimFile, ImageFile, saveout)
% Use to be 'FormatExperiment': Format data to be analyzed

saveInputData = false;
saveTriggers = false;
saveStimuli = false;
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
    saveout = true;
    SaveFile = StimFile;
elseif ischar(saveout)
    SaveFile = savout;
    saveout = true;
elseif islogical(saveout) && saveout == true
    SaveFile = StimFile;
end

%% Determine DataIn File
[p,fn,~] = fileparts(StimFile);
DataInFile = fullfile(p,strcat(fn,'_datain.bin')); %[strtok(StimFile, '.'),'_datain.bin'];

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

[~,~,ext] = fileparts(ImageFile);
switch ext
    case '.sbx'
        ImageInfo = parseSbxHeader(ImageFile);
    case '.imgs'
        ImageInfo = parseImgsHeader(ImageFile);
end

%% Load Data

% Load Data Files
load(StimFile, 'ImageInfo', 'StimuliNames', 'Triggers', 'Stimuli', 'DAQChannels', 'samplingFrequency', 'TrialInfo', 'Date', '-mat');
% info = fixTimeStamps(InfoFile); %for data collected w/ Scanbox v1.2
% info.event_id(info.event_id~=1) = []; %for data collected w/ Scanbox v1.2
load(InfoFile, 'info');
if ischar(DataInFile)
    DataInFID = fopen(DataInFile, 'r');
    DataIn = fread(DataInFID);
    fclose(DataInFID);
end

% Acquire Information from Image Files
if isempty(ImageInfo)
    error('user quit during parseSbxHeader. sbxPostProcess requires image files to be processed first');
end

% Acquire Information from Stim Files
nStimuliProtocols = length(StimuliNames);
nOutputChannels = size(Triggers, 2);
nInputChannels = length(DAQChannels) - nOutputChannels;

% Format DataIn
if ischar(DataInFile) && nInputChannels > 1 % if more than one channel is saved
    DataIn = permute(reshape(DataIn, nInputChannels, length(DataIn/nInputChannels)),[2,1]);
end


%% Creat Table With Row For Each Frame
frameRate = ImageInfo.FrameRate(1);
scansPerFrame = ceil(samplingFrequency/frameRate);
InputNames = DAQChannels(end-nInputChannels+1:end);
if saveout
    save(StimFile,'scansPerFrame','frameRate','InputNames','-append');
    VariablesToSave = {'AnalysisInfo', 'scansPerFrame', 'InputNames', 'frames'};
end

if any(strcmp(StimuliNames, 'Motor'))
    extra = rem(numel(info.event_id),4);
    if extra ~= 0 % imaging session started in middle of a trial (received odd # of triggers) so remove starting triggers
        info.event_id(1:extra) = [];
        info.frame(1:extra) = [];
        info.line(1:extra) = [];
    end
    nTrials = size(TrialInfo, 1);
    nTrials_imaging = length(info.event_id)/4;
    StimulusIndex = find(strcmp(StimuliNames, 'Scanbox')); % column in 'Stimuli' that matters for frames.Stimulus
elseif any(strcmp(StimuliNames, 'Piezo'))
    extra = rem(numel(info.event_id),2);
    if extra ~= 0 % imaging session started in middle of a trial (received odd # of triggers) so remove starting triggers
        info.event_id(1:extra) = [];
        info.frame(1:extra) = [];
        info.line(1:extra) = [];
    end
    nTrials = size(TrialInfo, 1);
    nTrials_imaging = length(info.event_id)/2;
    StimulusIndex = find(strcmp(StimuliNames, 'Piezo')); % column in 'Stimuli' that matters for frames.Stimulus
end
nTrialsMissed = nTrials - nTrials_imaging;
IndexStim = 1:nTrials;
IndexImaging = zeros(nTrials,1);
IndexImaging(IndexStim(nTrialsMissed+1:end)) = 1:nTrials_imaging;

AnalysisInfo = table(zeros(nTrials,1),zeros(nTrials,1),zeros(nTrials,1),cell(nTrials,1),zeros(nTrials,1),zeros(nTrials,1),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2), zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),...
    'VariableNames', {'StimID', 'TrialIndex', 'ImgIndex', 'ImgFilename', 'nFrames', 'nScans', 'ExpStimFrames', 'ExpStimScans', 'StimFrameLines', 'EventIDs', 'TrialStimFrames', 'TrialStimScans', 'ExpFrames', 'ExpScans'});

frames = struct('Stimulus', nan(ImageInfo.Frames(1),1), 'Trial', nan(ImageInfo.Frames(1),1));

if ismember(InputNames, 'RunSpeed')
    AnalysisInfo = [AnalysisInfo, table(cell(nTrials,1),'VariableNames',{'meanRunningSpeed'})];
    frames.RunningSpeed = nan(ImageInfo.Frames(1),1);
end
if saveInputData   
    AnalysisInfo = [AnalysisInfo, table(cell(nTrials,1),'VariableNames',{'InputData'})];
    VariablesToSave = [VariablesToSave, {'InputNames'}];
end
if saveTriggers
    AnalysisInfo = [AnalysisInfo, table(cell(nTrials,1),'VariableNames',{'Triggers'})];
    VariablesToSave = [VariablesToSave, {'DAQChannels'}];
end
if saveStimuli
    AnalysisInfo = [AnalysisInfo, table(cell(nTrials,1),'VariableNames',{'Stimuli'})];
    VariablesToSave = [VariablesToSave, {'StimuliNames'}];
end
if ischar(numBaselineFrames) || numBaselineFrames > 0
    AnalysisInfo = [AnalysisInfo, table(zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),zeros(nTrials,2),'VariableNames',{'ExpBaselineFrames','ExpBaselineScans','TrialBaselineFrames','TrialBaselineScans'})];
end

for t = 1:nTrials
    AnalysisInfo.StimID(t) = TrialInfo.TrialIndex(t); % ID # of stimulus presented
    AnalysisInfo.TrialIndex(t) = t; % index of trial presented relative to 'TrialInfo'
    AnalysisInfo.ImgIndex(t) = IndexImaging(t); % index of trial presented relative to what the images captured
    AnalysisInfo.ImgFilename{t} = ImageInfo.FullFilename(1); % filename of imaging file
    AnalysisInfo.nScans(t) = TrialInfo.NumScans(t); % number of scans in trial
    if t ~= 1 % determine index of scan that occurred last in previous trial
        lastscan = sum(TrialInfo.NumScans(1:t-1));
    else % first trial
        lastscan = 0;
    end
    AnalysisInfo.ExpScans(t,:) = [lastscan + 1, lastscan + AnalysisInfo.nScans(t)]; % indices of scans for current trial relative to entire experiment
    if t > nTrialsMissed
        if any(strcmp(StimuliNames, 'Motor'))
            AnalysisInfo.ExpStimFrames(t,:) = info.frame(4*(IndexImaging(t)-1)+2:4*IndexImaging(t)-1)'; % indices of first and last frame in which the stimulus occurred relative to entire experiment
            AnalysisInfo.StimFrameLines(t,:) = info.line(4*(IndexImaging(t)-1)+2:4*IndexImaging(t)-1)'; % indices of first and last line in which the stimulus occurred relative to each frame
            AnalysisInfo.EventIDs(t,:) = info.event_id(4*(IndexImaging(t)-1)+2:4*IndexImaging(t)-1)'; % ID # of TTL port in ScanBox that first and last trigger came in on
        elseif any(strcmp(StimuliNames, 'Piezo'))
            AnalysisInfo.ExpStimFrames(t,:) = info.frame(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % indices of first and last frame in which the stimulus occurred relative to entire experiment
            AnalysisInfo.StimFrameLines(t,:) = info.line(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % indices of first and last line in which the stimulus occurred relative to each frame
            AnalysisInfo.EventIDs(t,:) = info.event_id(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % ID # of TTL port in ScanBox that first and last trigger came in on
        end
        frames.Stimulus(AnalysisInfo.ExpStimFrames(t,1):AnalysisInfo.ExpStimFrames(t,2)) = AnalysisInfo.StimID(t);
        AnalysisInfo.TrialStimScans(t,:) = [find(Stimuli(:,StimulusIndex,AnalysisInfo.StimID(t)),1),find(Stimuli(:,StimulusIndex,AnalysisInfo.StimID(t)),1,'last')]; % indices of first and last scan in which the stimulus occurred relative to current trial
        AnalysisInfo.TrialStimFrames(t,:) = [floor(AnalysisInfo.TrialStimScans(t,1)/scansPerFrame),ceil(AnalysisInfo.TrialStimScans(t,2)/scansPerFrame)]; % indices of first and last frame in which the stimulus occurred relative to current trial
        AnalysisInfo.ExpFrames(t,:) = [AnalysisInfo.ExpStimFrames(t,1) - AnalysisInfo.TrialStimFrames(t,1) + 1, min(AnalysisInfo.ExpStimFrames(t,1) - AnalysisInfo.TrialStimFrames(t,1) + ceil(TrialInfo.TrialDuration(t)*frameRate),ImageInfo.Frames(1))]; % indices of frames for current trial relative to entire experiment
        frames.Trial(AnalysisInfo.ExpFrames(t,1):AnalysisInfo.ExpFrames(t,2)) = t;
        AnalysisInfo.nFrames(t) = AnalysisInfo.ExpFrames(t,2) - AnalysisInfo.ExpFrames(t,1) + 1;
        AnalysisInfo.ExpStimScans(t,:) = lastscan + AnalysisInfo.TrialStimScans(t,:);
    end
    if saveInputData
        AnalysisInfo.InputData{t} = DataIn(AnalysisInfo.ExpScans(t,1):min(end,AnalysisInfo.ExpScans(t,2)),:);
    end
    if ismember(InputNames, 'RunSpeed')
        AnalysisInfo.meanRunningSpeed{t} = zeros(AnalysisInfo.nFrames(t),1);
        index = find(strcmp(InputNames, 'RunSpeed'));
        if isfield(AnalysisInfo,'InputData')
            data = AnalysisInfo.InputData{t}(index,:);
        else
            data = DataIn(AnalysisInfo.ExpScans(t,1):min(end,AnalysisInfo.ExpScans(t,2)),index);
        end
        if ~isempty(data)
            data = calcRunningSpeed(data, samplingFrequency, 10, Date);
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
    if saveTriggers
        AnalysisInf.Triggers{t} = Triggers(:,:,AnalysisInfo.StimID(t));
    end
    if saveStimuli
        AnalysisInfo.Stimuli{t} = Stimuli(:,:,AnalysisInfo.StimID(t));
    end
    if ischar(numBaselineFrames) || numBaselineFrames > 0
        if ischar(numBaselineFrames)
            AnalysisInfo.ExpBaselineFrames(t,:) = [AnalysisInfo.ExpFrames(t,1), AnalysisInfo.ExpStimFrames(t,1)-1];
            AnalysisInfo.ExpBaselineScans(t,:) = [AnalysisInfo.ExpScans(t,1), AnalysisInfo.ExpStimScans(t,1) - 1];
            AnalysisInfo.TrialBaselineFrames(t,:) = [1, AnalysisInfo.TrialStimFrames(t,1) - 1];
            AnalysisInfo.TrialBaselineScans(t,:) = [1, AnalysisInfo.TrialStimScans(t,1) - 1];
        else
            AnalysisInfo.ExpBaselineFrames(t,:) = [max(1,AnalysisInfo.ExpStimFrames(t,1) - numBaselineFrames), AnalysisInfo.ExpStimFrames(t,1)-1];
            AnalysisInfo.ExpBaselineScans(t,:) = [max(1,AnalysisInfo.ExpStimScans(t,1) - numBaselineFrames*scansPerFrame), AnalysisInfo.ExpStimScans(t,1) - 1];
            AnalysisInfo.TrialBaselineFrames(t,:) = [max(1,AnalysisInfo.TrialStimFrames(t,1) - numBaselineFrames), AnalysisInfo.TrialStimFrames(t,1) - 1];
            AnalysisInfo.TrialBaselineScans(t,:) = [max(1,AnalysisInfo.TrialStimScans(t,1) - numBaselineFrames*scansPerFrame), AnalysisInfo.TrialStimScans(t,1) - 1];
        end
    end
end

% frames = struct2table(frames);

if saveout
    save(SaveFile, VariablesToSave{:}, '-append');
end

