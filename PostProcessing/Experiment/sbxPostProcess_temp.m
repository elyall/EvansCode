function [AnalysisInfo,frames,InputNames,Config] = sbxPostProcess_temp(ExperimentFile, ImageFile, WTFile, varargin)

directory = cd;

saveOut = false;
saveFile = '';


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
    directory = p;
end

if ~exist('ImageFile', 'var') || isempty(ImageFile)
    [ImageFile,p] = uigetfile({'*.sbx;*.tif;*.imgs'}, 'Choose corresponding Image File to process', directory);
    if isnumeric(ImageFile) % no file selected
        return
    end
    ImageFile = fullfile(p,ImageFile);
end

if ~exist('WTFile', 'var')
    WTFile = [];
end

if saveOut && isempty(saveFile)
    saveFile = ExperimentFile;
end

warning('off', 'MATLAB:table:RowsAddedExistingVars');
fprintf('Post-processing experiment file: %s\n', ExperimentFile);


%% Load in experiment data

% Load experiment info
load(ExperimentFile, 'DAQChannels', 'Experiment', 'TrialInfo', '-mat');
OutputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'O_')));
InputNames = DAQChannels(~cellfun(@isempty,strfind(DAQChannels, 'I_')));
nInputChannels = numel(InputNames);

% Load binary DAQ data
[p,fn,~] = fileparts(ExperimentFile);
DataInFile = fullfile(p, strcat(fn,'.bin')); %[strtok(StimFile, '.'),'_datain.bin'];
DataInFID = fopen(DataInFile, 'r');
DataIn = fread(DataInFID, inf, Experiment.saving.dataPrecision);
fclose(DataInFID);
DataIn = reshape(DataIn, nInputChannels, numel(DataIn)/nInputChannels)';


%% Load in neuronal imaging data
Config = load2PConfig(ImageFile);
[p,fn,~] = fileparts(ImageFile);
InfoFile = fullfile(p, strcat(fn,'.mat'));
load(InfoFile, 'info', '-mat');
info.frame = info.frame + 1;    % change 0 indexing to 1 indexing
info.line = info.line + 1;      % change 0 indexing to 1 indexing


%% Initialize outputs

AnalysisInfo = table(zeros(1,1),zeros(1,1),  cell(1,1),    cell(1,1),   zeros(1,1),zeros(1,2),zeros(1,2),   zeros(1,1), cell(1,1),   zeros(1,1),zeros(1,2),zeros(1,2),     zeros(1,2),       zeros(1,1),...
    'VariableNames', {'StimID','TrialIndex','ExpFilename','DataFilename','nScans','ExpScans','ExpStimScans','ImgIndex','ImgFilename','nFrames','ExpFrames','ExpStimFrames','StimFrameLines','onlineRunSpeed'});

frames = struct('Stimulus',nan(Config.Frames,1), 'Trial',nan(Config.Frames,1));

if isfield(TrialInfo, 'numRandomScansPost')
    toggleRandomTime = true;
    AnalysisInfo = [AnalysisInfo, table(zeros(1,1),'VariableNames',{'numRandomScansPost'})];
else
    toggleRandomTime = false;
end


%% Gather trial data

% Determine trial parameters
numTrials = numel(TrialInfo.StimID);
temp = Experiment.Triggers(:,strcmp(OutputNames,'O_2PTrigger'),1);
numTrigsPerTrial = nnz((temp-[0;temp(1:end-1)])>0);
numTrials2 = numel(info.frame(info.event_id==1))/numTrigsPerTrial;
if numTrials ~= numTrials2
    warning('Number of trials found in the experiment file do not match number of trials recorded by imaging computer!');
end

% Prepare trial assocation data
FrameCounter = DataIn(:,strcmp(InputNames,'I_FrameCounter'));
FrameCounter = FrameCounter - FrameCounter([1,1:end-1]);
FrameCounter = FrameCounter>0;

% Index through trials gathering necessary information
lastScan = 0;
for tindex = 1:numTrials
    
    % Assign trial info
    AnalysisInfo.StimID(tindex,1) = TrialInfo.StimID(tindex);             % ID # of stimulus presented
    AnalysisInfo.TrialIndex(tindex) = tindex;                           % index of trial from current '.exp' file
    AnalysisInfo.ExpFilename{tindex} = ExperimentFile;                  % name of '.exp' file current trial is found in
    AnalysisInfo.DataFilename{tindex} = DataInFile;                     % name of '.bin' file current input data is found in
    AnalysisInfo.onlineRunSpeed(tindex) = TrialInfo.RunSpeed(tindex);   % running speed calculated online
    if toggleRandomTime
        AnalysisInfo.numRandomScansPost(tindex) = TrialInfo.numRandomScansPost(tindex); % number of blank scans presented after current trial
    end
    
    % Determine scan info
    AnalysisInfo.nScans(tindex) = numel(Experiment.Stimulus);                                                           % number of scans in trial
    AnalysisInfo.ExpScans(tindex,:) = [lastScan+1, lastScan+AnalysisInfo.nScans(tindex)];                               % scans correlated to current trial
    AnalysisInfo.ExpStimScans(tindex,:) = lastScan + [find(Experiment.Stimulus,1),find(Experiment.Stimulus,1,'last')];  % scans in DataIn that current trial corresponds to (ignoring Random Scans Post)
    if ~toggleRandomTime
        lastScan = lastScan + AnalysisInfo.nScans(tindex);                                                              % update scan counter
    else
        lastScan = lastScan + AnalysisInfo.nScans(tindex) + AnalysisInfo.numRandomScansPost(tindex);                    % update scan counter
    end
    
    % Assign imaging info
    AnalysisInfo.ImgFilename{tindex} = Config.FullFilename; % filename of imaging file
    AnalysisInfo.ImgIndex(tindex) = tindex; % index of trial presented relative to what the images captured
    if numTrigsPerTrial == 2
        indices = numTrigsPerTrial*tindex-1:numTrigsPerTrial*tindex;
    elseif numTrigsPerTrial == 4
        indices = numTrigsPerTrial*tindex-2:numTrigsPerTrial*tindex-1;
    end
    if AnalysisInfo.StimID(tindex) ~= 0
        AnalysisInfo.ExpStimFrames(tindex,:) = info.frame(indices)'; % indices of frame in which trigger was received
        AnalysisInfo.StimFrameLines(tindex,:) = info.line(indices)'; % indices of line in which trigger was received
    else
        data = DataIn(AnalysisInfo.ExpStimScans(tindex-1,2):AnalysisInfo.ExpStimScans(tindex-1,2)+numel(Experiment.Stimulus),strcmp(InputNames,'I_FrameCounter'));
        data = data - [0;data(1:end-1)];
        data = data(2:end)>0;
        AnalysisInfo.ExpStimFrames(tindex,:) = AnalysisInfo.ExpFrames(tindex-1,2)+[sum(data(1:find(Experiment.Stimulus,1))), sum(data(1:find(Experiment.Stimulus,1,'last')))];
        AnalysisInfo.StimFrameLines(tindex,:) = nan(1,2);
        info.frame = [info.frame(1:indices(1)-1);AnalysisInfo.ExpStimFrames(tindex,:)';info.frame(indices(1):end)];
        info.line = [info.line(1:indices(1)-1);nan(2,1);info.line(indices(1):end)];
    end
    frames.Stimulus(AnalysisInfo.ExpStimFrames(tindex,1):AnalysisInfo.ExpStimFrames(tindex,2)) = AnalysisInfo.StimID(tindex); % define when stimulus was presented
    
    % Assign frames to individual trials
    if tindex == 1
        frameOffset = AnalysisInfo.ExpStimFrames(tindex,1) - sum(FrameCounter(1:AnalysisInfo.ExpStimScans(tindex,1))) - 1; %trig frame minus number of frames counted, minus 1 to account for frame that started before the experiment but ended within the experiment
    end
    AnalysisInfo.ExpFrames(tindex,:) = frameOffset + [sum(FrameCounter(1:AnalysisInfo.ExpScans(tindex,1))), sum(FrameCounter(1:AnalysisInfo.ExpScans(tindex,2)))]; % indices of frames for current trial relative to entire experiment
    AnalysisInfo.nFrames(tindex) = AnalysisInfo.ExpFrames(tindex,2) - AnalysisInfo.ExpFrames(tindex,1) + 1;
    if ~isnan(frames.Trial(AnalysisInfo.ExpFrames(tindex,1))) %first frame already associated with previous trial
        frames.Trial(AnalysisInfo.ExpFrames(tindex,1)+1:AnalysisInfo.ExpFrames(tindex,2)) = tindex; %do not overwrite first frame's association (because stimuli occur at end of trial)
    else
        frames.Trial(AnalysisInfo.ExpFrames(tindex,1):AnalysisInfo.ExpFrames(tindex,2)) = tindex;
    end
    
end %trials


%% Determine alignment of 2P frames to scans
data = DataIn(:,strcmp(InputNames,'I_FrameCounter'));
data = data - data([1,1:end-1]);
frames.StartScan = nan(Config.Frames,1);
if data(1)==1
    data = find(data>0);
    frames.StartScan(AnalysisInfo.ExpFrames(1,1):AnalysisInfo.ExpFrames(1,1)+numel(data)-1) = data;
else
    data = find(data>0);
    frames.StartScan(AnalysisInfo.ExpFrames(1,1)+1:AnalysisInfo.ExpFrames(1,1)+numel(data)) = data;
end
% currentFrame = AnalysisInfo.ExpFrames(1,1);
% for findex=1:numel(ScanIndex)
%     frames.ScanIndex(currentFrame) = ScanIndex(findex);
%     currentFrame = currentFrame + 1;
% end


%% Compute running speed

% Compute running speed at acquired sampling frequency
data = DataIn(:,strcmp(InputNames, 'I_RunWheelB'));
dataB = DataIn(:,strcmp(InputNames, 'I_RunWheelA'));
data = data - [0;data(1:end-1)];
data = data>0;
data(all([data,dataB],2))=-1; %set backwards steps to be backwards
data = calcRunningSpeed(data, Experiment.params.samplingFrequency, 1, date);

% Set running speed per frame to be mean during that frame
frames.RunningSpeed = nan(Config.Frames,1);
frames.RunningSpeed(AnalysisInfo.ExpFrames(1,1)) = mean(data(1:frames.StartScan(AnalysisInfo.ExpFrames(1,1)+1)-1));
for frameNumber = AnalysisInfo.ExpFrames(1,1)+1:find(~isnan(frames.StartScan),1,'last')-1
    frames.RunningSpeed(frameNumber) = mean(data(frames.StartScan(frameNumber):frames.StartScan(frameNumber+1)-1));
end


%% Determine alignment of Whisker Tracking Frames
if any(ismember(InputNames, 'I_WhiskerTracker'))
    if exist(WTFile, 'file')
        % Determine number of frames recorded
        temp = VideoReader(WTFile);
        numFrames = temp.Duration*temp.FrameRate;
        
        % Determine number of frames captured
        data = DataIn(:,strcmp(InputNames, 'I_WhiskerTracker'));
        data = data - [0;data(1:end-1)];
        data = find(data>0);
        numFrames2 = numel(data);
        
        if numFrames < numFrames2
            warning('Fuck: more WT frames registered than actually recorded');
        end
        
        % Assign frames to specific trigger
        AnalysisInfo = [AnalysisInfo, table(cell(numTrials,1),'VariableNames',{'Scan2WTFrame'})];
        frameCount = 0;
        for tindex = 1:numTrials
            indices = find(data>=AnalysisInfo.ExpScans(tindex,1) & data<=AnalysisInfo.ExpScans(tindex,2));
            AnalysisInfo.Scan2WTFrame{tindex} = [data(indices),(frameCount+1:frameCount+numel(indices))'];
            frameCount = frameCount + numel(indices);
        end
    else
        warning('Whisker Images file not found, will not include WT data');
    end
end


%% Save output
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'AnalysisInfo', 'frames', 'InputNames', '-mat', '-v7.3');
    else
        save(saveFile, 'AnalysisInfo', 'frames', 'InputNames', '-mat', '-append');
    end
    fprintf('AnalysisInfo saved to: %s\n', saveFile);
end

