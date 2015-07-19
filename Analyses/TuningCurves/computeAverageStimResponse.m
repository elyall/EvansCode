function [StimResponse, AvgFrames, AvgDFoFFrames] = computeAverageStimResponse(ImageFile, ExperimentFile, minRunSpeed, MotionCorrect)


%% Check input arguments
narginchk(0, 4)
if ~exist('ImageFile','var') || isempty(ImageFile) % Prompt for file selection
    directory = CanalSettings('DataDirectory');
    [ImageFile,p] = uigetfile({'*.imgs;*.sbx;*.tif'},'Select images:',directory);
    if isnumeric(ImageFile)
        return
    end
    ImageFile = fullfile(p, ImageFile);
end

if ~exist('ExperimentFile','var') || isempty(ExperimentFile)
    directory = CanalSettings('ExperimentDirectory');
    [ExperimentFile, p] = uigetfile({'*.mat'},'Choose Experiment file',directory);
    if isnumeric(ExperimentFile)
        return
    end
    ExperimentFile = fullfile(p,ExperimentFile);
    load(ExperimentFile, 'AnalysisInfo', '-mat');
elseif ischar(ExperimentFile)
    load(ExperimentFile, 'AnalysisInfo', '-mat');
elseif isstruct(ExperimentFile)
    AnalysisInfo = ExperimentFile;
    clear ExperimentFile
end

if ~exist('minRunSpeed', 'var') || isempty(minRunSpeed)
    minRunSpeed = 100;
end

if ~exist('MotionCorrect', 'var') || isempty(MotionCorrect)
    vars = whos(matfile(ExperimentFile));
    if any(strcmp(vars, 'MCdata'))
        MotionCorrect = true;
    else
        MotionCorrect = false;
    end
end

%% Motion correction
if MotionCorrect
    load(ExperimentFile, 'MCdata', '-mat');
end

%% Load in image config
fprintf('\nCalculating average stimuli for: %s\n', ImageFile);
[~, Config] = load2P(ImageFile,'Type','Direct','Frames',2);

%% Determine frames to load for each stimulus
StimIDs = unique(AnalysisInfo.StimID);
nStims = numel(StimIDs);
StimResponse = struct();
AvgFrames = cell(nStims,1);
AvgDFoFFrames = cell(nStims,1);
AvgEvokedDFoF = zeros([Config.size(1:end-2), nStims]);
SEEvokedDFoF =  zeros([Config.size(1:end-2), nStims]);
EvokedDFoFPValue = zeros([Config.size(1:end-2), nStims-1]); %significantly different response relative to control trials
Excited = zeros([Config.size(1:end-2), nStims-1]); %increased activity relative to control trials
for sindex = 1:nStims
    StimResponse(sindex).StimID = StimIDs(sindex); % stimulus ID
    StimResponse(sindex).TrialIDs = find(AnalysisInfo.StimID==StimIDs(sindex)); % trial IDs
    StimResponse(sindex).nTrials = numel(StimResponse(sindex).TrialIDs);
    StimResponse(sindex).trialStimFrames = AnalysisInfo.TrialStimFrames(StimResponse(sindex).TrialIDs,:); % stimulus frames relative to trial
    StimResponse(sindex).nStimFrames = diff(StimResponse(sindex).trialStimFrames,1,2)+1;
    StimResponse(sindex).trialBaselineFrames = AnalysisInfo.TrialBaselineFrames(StimResponse(sindex).TrialIDs,:); % baseline frames relative to trial
    StimResponse(sindex).expStimFrames = AnalysisInfo.ExpStimFrames(StimResponse(sindex).TrialIDs,:); % stimulus frames relative to experiment
    StimResponse(sindex).expFrames = AnalysisInfo.ExpFrames(StimResponse(sindex).TrialIDs,:); % frames for current trial
    StimResponse(sindex).nFrames = AnalysisInfo.nFrames(StimResponse(sindex).TrialIDs,:); % number of frames in each trial
    StimResponse(sindex).imaged = logical(AnalysisInfo.ImgIndex(StimResponse(sindex).TrialIDs)); % whether or not the trial was captured by the pmts
    StimResponse(sindex).imaged(StimResponse(sindex).nFrames~=mode(StimResponse(sindex).nFrames)) = false;
    StimResponse(sindex).GoodTrials = StimResponse(sindex).imaged;
    if isfield(AnalysisInfo, 'meanRunningSpeed')
        StimResponse(sindex).RunSpeed = AnalysisInfo.meanRunningSpeed(StimResponse(sindex).TrialIDs); % running speed for each trial
        temp = [StimResponse(sindex).RunSpeed{StimResponse(sindex).imaged}];
        StimResponse(sindex).running = StimResponse(sindex).GoodTrials;
        StimResponse(sindex).running(StimResponse(sindex).imaged) = mean(temp, 1) >= minRunSpeed;
        StimResponse(sindex).GoodTrials = StimResponse(sindex).running; % number of imaged trials in which mouse was running for current stimulus
    end
    StimResponse(sindex).nGoodTrials = sum(StimResponse(sindex).GoodTrials);
        
    % Compute Average Stimulus
    nFrames = mode(StimResponse(sindex).nFrames);
    nStimFrames = mode(StimResponse(sindex).nStimFrames);
    AvgFrames{sindex} = zeros([Config.size(1:end-1), nFrames]);
    AvgDFoFFrames{sindex} = zeros([Config.size(1:end-2), nFrames]);
    StimResponse(sindex).trialEvokedDFoF = zeros([Config.size(1:end-2), StimResponse(sindex).nTrials]);
    [Frames, ~, loadObj] = load2P(ImageFile, 'Type', 'MemMap', 'Double'); %memmap
    if strcmp(loadObj.ext{1}, '.sbx') && strcmp(loadObj.Type, 'MemMap') %memap
        Frames = intmax(Config.Precision) - Frames; %memmap
    end %memmap
    % Load in one trial at at time
    for tindex = 1:StimResponse(sindex).nTrials
        if StimResponse(sindex).imaged(tindex)
            frames = double(Frames(:,:,:,:,StimResponse(sindex).expFrames(tindex,1):StimResponse(sindex).expFrames(tindex,2))); %memmap
%             frames = load2P(ImageFile, 'Frames', StimResponse(sindex).expFrames(tindex,1):StimResponse(sindex).expFrames(tindex,2), 'Double'); %direct
            if MotionCorrect
                frames = applyMotionCorrection(frames, MCdata, StimResponse(sindex).expFrames(tindex,1):StimResponse(sindex).expFrames(tindex,2));
            end
            baselineFrame = mean(frames(:,:,:,1,StimResponse(sindex).trialBaselineFrames(tindex,1):StimResponse(sindex).trialBaselineFrames(tindex,2)),5);
            baselineStim = repmat(baselineFrame,1,1,Config.size(3),nStimFrames);
            StimResponse(sindex).trialEvokedDFoF(:,:,:,tindex) = mean((permute(frames(:,:,:,1,StimResponse(sindex).trialStimFrames(tindex,1):StimResponse(sindex).trialStimFrames(tindex,2)),[1,2,3,5,4])-baselineStim)./baselineStim,4);
            if StimResponse(sindex).GoodTrials(tindex)
                AvgFrames{sindex} = AvgFrames{sindex} + frames/StimResponse(sindex).nGoodTrials;
                baselineAvg = repmat(baselineFrame,1,1,Config.size(3),nFrames);
                AvgDFoFFrames{sindex} = AvgDFoFFrames{sindex} + ((permute(frames(:,:,:,1,:),[1,2,3,5,4])-baselineAvg)./baselineAvg)/StimResponse(sindex).nGoodTrials;
            end
        end
    end  
    AvgEvokedDFoF(:,:,:,sindex) = mean(StimResponse(sindex).trialEvokedDFoF(:,:,:,StimResponse(sindex).GoodTrials), 4);
    SEEvokedDFoF(:,:,:,sindex) = std(StimResponse(sindex).trialEvokedDFoF(:,:,:,StimResponse(sindex).GoodTrials),[],4)/sqrt(StimResponse(sindex).nGoodTrials);
    if sindex ~= 1
        [~, EvokedDFoFPValue(:,:,:,sindex-1)] = ttest2(...
            StimResponse(1).trialEvokedDFoF(:,:,:,StimResponse(1).GoodTrials),....
            StimResponse(sindex).trialEvokedDFoF(:,:,:,StimResponse(sindex).GoodTrials),...
            'dim', 4);
        Excited(:,:,:,sindex-1) = AvgEvokedDFoF(:,:,:,sindex) >= AvgEvokedDFoF(:,:,:,1);
    end
    fprintf('\tstim %d of %d: calculated response (%d trials)\n',sindex,nStims,StimResponse(sindex).nGoodTrials);
end

%% Save to file
if exist('ExperimentFile', 'var') && ischar(ExperimentFile)
    save(ExperimentFile, 'AvgEvokedDFoF', 'SEEvokedDFoF', 'StimResponse', 'EvokedDFoFPValue', 'AvgFrames', 'AvgDFoFFrames', '-append');
    fprintf('Saved average stimuli to: %s\n', ExperimentFile);
end