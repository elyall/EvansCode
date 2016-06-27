function [AnalysisInfo,frames,InputNames,Config] = sbxPostProcess2(ExperimentFiles, ImageFiles, varargin)
% Use to be 'FormatExperiment': Format data to be analyzed

directory = cd;

saveOut = false;
saveFile = '';
saveInputData = false;

numTrigsPerTrial = 2;

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

if ~exist('ExperimentFiles', 'var') || isempty(ExperimentFiles) % nothing input
    [ExperimentFiles,p] = uigetfile({'*.exp;*.mat'},'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ExperimentFiles) % no file selected
        return
    elseif iscellstr(ExperimentFiles) % multiple files selected
        ExperimentFiles = fullfile(p, ExperimentFiles);
    elseif ischar(ExperimentFiles) % single file selected
        ExperimentFiles = {fullfile(p, ExperimentFiles)};
    end
elseif ischar(ExperimentFiles)
    if isdir(ExperimentFiles) % directory input
        p = ExperimentFiles;
        ExperimentFiles = dir(fullfile(p, '*.exp'));
        ExperimentFiles = fullfile(p, {ExperimentFiles.name});
    else % single file input
        ExperimentFiles = {ExperimentFiles};
    end
end

if ~exist('ImageFiles', 'var') || isempty(ImageFiles)
    [ImageFiles,p] = uigetfile({'*.sbx;*.tif;*.imgs'}, 'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ImageFiles) % no file selected
        return
    elseif iscellstr(ImageFiles) % multiple files selected
        ImageFiles = fullfile(p, ImageFiles);
    elseif ischar(ImageFiles) % single file selected
        ImageFiles = {fullfile(p, ImageFiles)};
    end
elseif ischar(ImageFiles)
    if isdir(ImageFiles) % directory input
        p = ImageFiles;
        ImageFiles = dir(p);
        ImageFiles = fullfile(p, {ImageFiles(~cellfun(@isempty, regexpi({ImageFiles.name}, '.*(sbx|tif)'))).name});
    else % single file input
        ImageFiles = {ImageFiles};
    end
end

if saveOut && isempty(saveFile)
    saveFile = ExperimentFiles{1};
end

numFiles = numel(ExperimentFiles);


%% Determine number of frames in experiment
Config = load2PConfig(ImageFiles);
totalFrames = sum([Config(:).Frames]);


%% Determine input data available
load(ExperimentFiles{1}, 'DAQChannels', '-mat');
OutputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'O_')));
InputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
nInputChannels = numel(InputNames);


%% Initialize trial table and series variable
if saveOut
    VariablesToSave = {'AnalysisInfo', 'frames', 'InputNames'};
end

AnalysisInfo = table(zeros(1,1),zeros(1,1),zeros(1,1),cell(1,1),cell(1,1),cell(1,1),zeros(1,1),zeros(1,1),zeros(1,2),zeros(1,2),zeros(1,2),zeros(1,2),zeros(1,2),zeros(1,2),zeros(1,2),zeros(1,1),...
    'VariableNames', {'StimID', 'TrialIndex', 'ImgIndex', 'ImgFilename', 'ExpFilename', 'DataFilename', 'nFrames', 'nScans', 'ExpStimFrames', 'ExpStimScans', 'StimFrameLines', 'TrialStimFrames', 'TrialStimScans', 'ExpFrames', 'ExpScans', 'scansPerFrame'});

frames = struct('Stimulus', nan(totalFrames,1), 'Trial', nan(totalFrames,1));

if any(ismember(InputNames, 'I_RunWheelA'))
    AnalysisInfo = [AnalysisInfo, table(cell(1,1),'VariableNames',{'meanRunningSpeed'})];
    frames.RunningSpeed = nan(totalFrames,1);
end

if saveInputData
    AnalysisInfo = [AnalysisInfo, table(cell(1,1),'VariableNames',{'InputData'})];
end

warning('off', 'MATLAB:table:RowsAddedExistingVars');
trialCounter = 0;
fprintf('Post-processing %d experiment file(s)...\n', numFiles);
fprintf('\t%s\n', ExperimentFiles{:});
for findex = 1:numFiles
    
    % Load Experiment Info
    load(ExperimentFiles{findex}, 'Experiment', 'TrialInfo', '-mat');
    [p,fn,~] = fileparts(ExperimentFiles{findex});
    DataInFile = fullfile(p, strcat(fn,'.bin')); %[strtok(StimFile, '.'),'_datain.bin'];
    scansPerFrame = Experiment.params.samplingFrequency/Config(findex).FrameRate;

    % Format DataIn
    if exist(DataInFile, 'file')
        DataInFID = fopen(DataInFile, 'r');
        DataIn = fread(DataInFID, inf, Experiment.saving.dataPrecision);
        fclose(DataInFID);
    else
        DataInFile = false;
    end
    if ischar(DataInFile) && nInputChannels > 1 % if more than one channel is saved
        DataIn = reshape(DataIn, nInputChannels, numel(DataIn)/nInputChannels)';
    end
    
    % Load in Image info
    [p,fn,~] = fileparts(ImageFiles{findex});
    InfoFile = fullfile(p, strcat(fn,'.mat')); %[strtok(StimFile, '.'),'_datain.bin'];
    load(InfoFile, 'info');
    info = fixTimeStamps(InfoFile); %for data collected w/ Scanbox v1.2
    info.frame(info.event_id~=1) = []; %remove info from second input channel
    info.line(info.event_id~=1) = [];
    info.event_id(info.event_id~=1) = [];
    info.frame = info.frame + 1; % change 0 indexing to 1 indexing
    
    % Determine trial parameters
    numTrials = numel(TrialInfo.StimID);
    temp = Experiment.Triggers(:,strcmp(OutputNames,'O_2PTrigger'),1);
    numTrigsPerTrial = nnz((temp-[0;temp(1:end-1)])>0);
    nTrials_imaging = numel(info.frame(info.event_id==1))/numTrigsPerTrial;
    nTrialsMissed = numTrials - nTrials_imaging;
    IndexStim = 1:numTrials;
    IndexImaging = zeros(numTrials,1);
    IndexImaging(IndexStim(nTrialsMissed+1:end)) = 1:nTrials_imaging;
    
    for tindex = 1:numTrials
        
        % Assign trial info
        AnalysisInfo.StimID(tindex+trialCounter,1) = TrialInfo.StimID(tindex); % ID # of stimulus presented
        AnalysisInfo.TrialIndex(tindex+trialCounter) = tindex; % index of trial presented relative to 'TrialInfo'
        AnalysisInfo.ImgIndex(tindex+trialCounter) = IndexImaging(tindex); % index of trial presented relative to what the images captured
        AnalysisInfo.ImgFilename{tindex+trialCounter} = Config(findex).FullFilename; % filename of imaging file
        AnalysisInfo.ExpFilename{tindex+trialCounter} = ExperimentFiles{findex};
        AnalysisInfo.DataFilename{tindex+trialCounter} = DataInFile;
        AnalysisInfo.scansPerFrame(tindex+trialCounter) = scansPerFrame;
        
        % Assign trial scans and frames
        AnalysisInfo.nScans(tindex+trialCounter) = numel(Experiment.Stimulus); % number of scans in trial
        if tindex ~= 1 % determine index of scan that occurred last in previous trial
            lastscan = sum(AnalysisInfo.nScans(1:tindex-1));
        else % first trial
            lastscan = 0;
        end
        AnalysisInfo.ExpScans(tindex+trialCounter,:) = [lastscan + 1, lastscan + AnalysisInfo.nScans(tindex)]; % indices of scans for current trial relative to entire experiment
        if tindex > nTrialsMissed
            AnalysisInfo.ExpStimFrames(tindex+trialCounter,:) = info.frame(4*(IndexImaging(tindex)-1)+2:4*IndexImaging(tindex)-1)'; % indices of first and last frame in which the stimulus occurred relative to entire experiment
            AnalysisInfo.StimFrameLines(tindex+trialCounter,:) = info.line(4*(IndexImaging(tindex)-1)+2:4*IndexImaging(tindex)-1)'; % indices of first and last line in which the stimulus occurred relative to each frame
            %AnalysisInfo.ExpStimFrames(t,:) = info.frame(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % indices of first and last frame in which the stimulus occurred relative to entire experiment
            %AnalysisInfo.StimFrameLines(t,:) = info.line(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % indices of first and last line in which the stimulus occurred relative to each frame
            %AnalysisInfo.EventIDs(t,:) = info.event_id(2*(IndexImaging(t)-1)+1:2*IndexImaging(t))'; % ID # of TTL port in ScanBox that first and last trigger came in on
            frames.Stimulus(AnalysisInfo.ExpStimFrames(tindex+trialCounter,1):AnalysisInfo.ExpStimFrames(tindex+trialCounter,2)) = AnalysisInfo.StimID(tindex+trialCounter);
            AnalysisInfo.TrialStimScans(tindex+trialCounter,:) = [find(Experiment.Stimulus,1),find(Experiment.Stimulus,1,'last')]; % indices of first and last scan in which the stimulus occurred relative to current trial
            AnalysisInfo.TrialStimFrames(tindex+trialCounter,:) = [floor(AnalysisInfo.TrialStimScans(tindex+trialCounter,1)/scansPerFrame),ceil(AnalysisInfo.TrialStimScans(tindex+trialCounter,2)/scansPerFrame)]; % indices of first and last frame in which the stimulus occurred relative to current trial
            AnalysisInfo.ExpFrames(tindex+trialCounter,:) = [AnalysisInfo.ExpStimFrames(tindex+trialCounter,1) - AnalysisInfo.TrialStimFrames(tindex+trialCounter,1) + 1, min(AnalysisInfo.ExpStimFrames(tindex+trialCounter,1) - AnalysisInfo.TrialStimFrames(tindex+trialCounter,1) + ceil(numel(Experiment.Stimulus)/scansPerFrame),Config(findex).Frames)]; % indices of frames for current trial relative to entire experiment
            frames.Trial(AnalysisInfo.ExpFrames(tindex+trialCounter,1):AnalysisInfo.ExpFrames(tindex+trialCounter,2)) = tindex+trialCounter;
            AnalysisInfo.nFrames(tindex+trialCounter) = AnalysisInfo.ExpFrames(tindex+trialCounter,2) - AnalysisInfo.ExpFrames(tindex+trialCounter,1) + 1;
            AnalysisInfo.ExpStimScans(tindex+trialCounter,:) = lastscan + AnalysisInfo.TrialStimScans(tindex,:);
        end
        
        % Save input data if requested
        if saveInputData
            AnalysisInfo.InputData{tindex+trialCounter} = DataIn(AnalysisInfo.ExpScans(tindex+trialCounter,1):min(end,AnalysisInfo.ExpScans(tindex+trialCounter,2)),:);
        end
        
        % Compute and save running speed information
        if any(ismember(InputNames, 'I_RunWheelA'))
            AnalysisInfo.meanRunningSpeed{tindex+trialCounter} = zeros(AnalysisInfo.nFrames(tindex+trialCounter),1);
            data = DataIn(AnalysisInfo.ExpScans(tindex+trialCounter,1):min(end,AnalysisInfo.ExpScans(tindex+trialCounter,2)),find(strcmp(InputNames, 'I_RunWheelB')));
            dataB = DataIn(AnalysisInfo.ExpScans(tindex+trialCounter,1):min(end,AnalysisInfo.ExpScans(tindex+trialCounter,2)),find(strcmp(InputNames, 'I_RunWheelA')));
            data = (data - [0; data(1:end-1)]);
            data(data<0) = 0;
            dataB(~data) = 0;
            data(logical(dataB)) = -1;
            if ~isempty(data)
                data = calcRunningSpeed(data, Experiment.params.samplingFrequency, 10, date);
                dataPointsPerFrame = round(scansPerFrame/10);
                for frameNumber = 1:AnalysisInfo.nFrames(tindex+trialCounter)
                    AnalysisInfo.meanRunningSpeed{tindex+trialCounter}(frameNumber) = mean(data(dataPointsPerFrame*(frameNumber-1)+1:min(end,dataPointsPerFrame*frameNumber)));
                    frames.RunningSpeed(AnalysisInfo.ExpFrames(tindex+trialCounter,1)+frameNumber-1) = AnalysisInfo.meanRunningSpeed{tindex+trialCounter}(frameNumber);
                end
            else
                for frameNumber = 1:AnalysisInfo.nFrames(tindex+trialCounter)
                    AnalysisInfo.meanRunningSpeed{tindex+trialCounter}(frameNumber) = NaN;
                    frames.RunningSpeed(AnalysisInfo.ExpFrames(tindex+trialCounter,1)+frameNumber-1) = AnalysisInfo.meanRunningSpeed{tindex+trialCounter}(frameNumber);
                end
            end
        end % running speed
        
    end % trials
    
    trialCounter = trialCounter + numTrials; % indexing trial counter
    
end % files


%% Save output
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, VariablesToSave{:}, '-mat', '-v7.3');
    else
        save(saveFile, VariablesToSave{:}, '-mat', '-append');
    end
    fprintf('AnalysisInfo saved to: %s\n', saveFile);
end

