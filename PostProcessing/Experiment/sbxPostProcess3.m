function [AnalysisInfo,frames,InputNames,Config] = sbxPostProcess3(ExperimentFile, ImageFile, WTFile, varargin)
% SBXPOSTPROCESS3 Process experiment data
%   sbxPostProcess3(ExperimentFile, ImageFile, WhiskerTrackingFile)
%   produces an AnalysisInfo table with the number of rows equal to the
%   number of trials recorded in the ExperimentFile. Each row has all the
%   relevant information for that trial.
%
% Update in version 3: computes frame wise information based on triggers
% captured per frame, rather than estimating off of the triggers sent per
% stimulus.

directory = cd;

saveOut = false;
saveFile = '';

downsample_runspeed = 1;

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
fprintf('Post-processing experiment file: %s', ExperimentFile);


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

AnalysisInfo = table(zeros(1,1),zeros(1,1),  cell(1,1),    cell(1,1),   zeros(1,1),zeros(1,2),zeros(1,2),   zeros(1,1), cell(1,1),   zeros(1,1),zeros(1,2),zeros(1,2),     zeros(1,2),...
    'VariableNames', {'StimID','TrialIndex','ExpFilename','DataFilename','nScans','ExpScans','ExpStimScans','ImgIndex','ImgFilename','nFrames','ExpFrames','ExpStimFrames','StimFrameLines'});

frames = struct('Stimulus',nan(Config.Frames,1), 'Trial',nan(Config.Frames,1));

if any(ismember(InputNames,'I_RunWheelB'))
    runspeed = true;
    AnalysisInfo = [AnalysisInfo, table(zeros(1,1),'VariableNames',{'onlineRunSpeed'})];
else
    runspeed = false;
end

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
if any(DataIn(:,strcmp(InputNames,'I_FrameCounter')))
    FrameCounter = DataIn(:,strcmp(InputNames,'I_FrameCounter'));
    FrameCounter = FrameCounter - FrameCounter([1,1:end-1]);
    FrameCounter = FrameCounter>0;
else
    warning('Frame counter failed to record any frames -> estimating exact frame cut offs');
    FrameCounter = false;
    scansPerFrame = Experiment.params.samplingFrequency/Config.FrameRate;
end

% Index through trials gathering necessary information
lastScan = 0;
for tindex = 1:numTrials
    
    % Assign trial info
    AnalysisInfo.StimID(tindex,1) = TrialInfo.StimID(tindex);           % ID # of stimulus presented
    AnalysisInfo.TrialIndex(tindex) = tindex;                           % index of trial from current '.exp' file
    AnalysisInfo.ExpFilename{tindex} = ExperimentFile;                  % name of '.exp' file current trial is found in
    AnalysisInfo.DataFilename{tindex} = DataInFile;                     % name of '.bin' file current input data is found in
    if runspeed
        AnalysisInfo.onlineRunSpeed(tindex,1) = TrialInfo.RunSpeed(tindex); % running speed calculated online
    end
    if toggleRandomTime
        AnalysisInfo.numRandomScansPost(tindex,1) = TrialInfo.numRandomScansPost(tindex); % number of blank scans presented after current trial
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
    AnalysisInfo.ExpStimFrames(tindex,:) = info.frame(indices)'; % indices of frame in which trigger was received
    AnalysisInfo.StimFrameLines(tindex,:) = info.line(indices)'; % indices of line in which trigger was received
    frames.Stimulus(AnalysisInfo.ExpStimFrames(tindex,1):AnalysisInfo.ExpStimFrames(tindex,2)) = AnalysisInfo.StimID(tindex); % define when stimulus was presented
    
    % Assign frames to individual trials
    if ~isequal(FrameCounter,false)
        if tindex == 1
            frameOffset = AnalysisInfo.ExpStimFrames(tindex,1) - sum(FrameCounter(1:AnalysisInfo.ExpStimScans(tindex,1))); %trig frame minus number of full frames recorded during the experiment
        end
        AnalysisInfo.ExpFrames(tindex,:) = frameOffset + [sum(FrameCounter(1:AnalysisInfo.ExpScans(tindex,1))), sum(FrameCounter(1:AnalysisInfo.ExpScans(tindex,2)))]; % indices of frames for current trial relative to entire experiment
        
    else
        N = round((AnalysisInfo.ExpScans(tindex,2)-AnalysisInfo.ExpStimScans(tindex,2))/scansPerFrame);
        if tindex == 1
            AnalysisInfo.ExpFrames(tindex,:) = [1,AnalysisInfo.ExpStimFrames(1,2)+N];
        else
            AnalysisInfo.ExpFrames(tindex,:) = [AnalysisInfo.ExpFrames(tindex-1,2)+1,AnalysisInfo.ExpStimFrames(tindex,2)+N];
        end
    end
    AnalysisInfo.nFrames(tindex) = AnalysisInfo.ExpFrames(tindex,2) - AnalysisInfo.ExpFrames(tindex,1) + 1;
    if ~isnan(frames.Trial(AnalysisInfo.ExpFrames(tindex,1))) %first frame already associated with previous trial
        frames.Trial(AnalysisInfo.ExpFrames(tindex,1)+1:AnalysisInfo.ExpFrames(tindex,2)) = tindex; %do not overwrite first frame's association (because stimuli occur at end of trial)
    else
        frames.Trial(AnalysisInfo.ExpFrames(tindex,1):AnalysisInfo.ExpFrames(tindex,2)) = tindex;
    end
    
end %trials


%% Determine alignment of 2P frames to scans
frames.StartScan = nan(Config.Frames,1);
if any(ismember(InputNames,'I_FrameCounter')) && any(DataIn(:,strcmp(InputNames,'I_FrameCounter')))
    data = DataIn(:,strcmp(InputNames,'I_FrameCounter'));
    data = data - data([1,1:end-1]);
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
end


%% Compute running speed
if runspeed
    frames.RunningSpeed = nan(Config.Frames,1);
    if any(DataIn(:,strcmp(InputNames,'I_RunWheelB'))) % cable was plugged in properly
        
        % Compute running speed at acquired sampling frequency
        data = DataIn(:,strcmp(InputNames, 'I_RunWheelB'));
        dataB = DataIn(:,strcmp(InputNames, 'I_RunWheelA'));
        data = (data - [0; data(1:end-1)]);
        data = data>0;
        data(all([data,dataB],2))=-1; %set backwards steps to be backwards
        data = calcRunningSpeed(data, Experiment.params.samplingFrequency, downsample_runspeed, date);
        
        % Determine mean running speed per frame
        if ~all(isnan(frames.StartScan)) % exact scans for each frame are known
            frames.RunningSpeed(AnalysisInfo.ExpFrames(1,1)) = mean(data(1:frames.StartScan(AnalysisInfo.ExpFrames(1,1)+1)-1));
            for frameNumber = AnalysisInfo.ExpFrames(1,1)+1:find(~isnan(frames.StartScan),1,'last')-1
                frames.RunningSpeed(frameNumber) = mean(data(frames.StartScan(frameNumber):frames.StartScan(frameNumber+1)-1));
            end
        else % FrameCounter didn't exist or was empty -> assume sampling to be constant
            for tindex = 1:numTrials
                temp = data(max(1,round(AnalysisInfo.ExpScans(tindex,1)/downsample_runspeed)):round(AnalysisInfo.ExpScans(tindex,2)/downsample_runspeed));
                N = floor(numel(temp)/AnalysisInfo.nFrames(tindex));
                temp = reshape(temp(1:N*AnalysisInfo.nFrames(tindex)),N,AnalysisInfo.nFrames(tindex));
                frames.RunningSpeed(AnalysisInfo.ExpFrames(tindex,1):AnalysisInfo.ExpFrames(tindex,2)) = mean(temp);
            end
        end
    end
end


%% Determine alignment of Whisker Tracking Frames
if any(ismember(InputNames, 'I_WhiskerTracker'))
    if exist(WTFile, 'file')
        % Determine number of frames recorded
        temp = VideoReader(WTFile);
        numFrames = temp.Duration*temp.FrameRate;
        
        % Determine number of frames registered
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
            indices = find(data>=AnalysisInfo.ExpScans(tindex,1) && data<=AnalysisInfo.ExpScans(tindex,2));
            AnalysisInfo.Scan2WTFrame{tindex} = [data(indices),(frameCount+1:frameCount+numel(indices))'];
            frameCount = frameCount + numel(indices);
        end
    else
        warning('Whisker Images file not found, will not include WT data');
    end
end
fprintf('\tComplete.\n');


%% Save output
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'AnalysisInfo', 'frames', 'InputNames', '-mat', '-v7.3');
    else
        save(saveFile, 'AnalysisInfo', 'frames', 'InputNames', '-mat', '-append');
    end
    fprintf('\tAnalysisInfo saved to: %s\n', saveFile);
end

