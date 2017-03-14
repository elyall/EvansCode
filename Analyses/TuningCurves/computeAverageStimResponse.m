function [AvgTrial, AvgTrialdFoF, StimResponse] = computeAverageStimResponse(ImageFiles, ExperimentFile, TrialIndex, MotionCorrect, varargin)

Depth = 1;
Channel = 1;
timeBefore = 2;
timeAfter = 4;

saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
relativeID(1) = 1;
while relativeID(1)<=length(varargin)
    try
        switch varargin{relativeID(1)}
            case {'Depth','depth'}
                Depth = varargin{relativeID(1)+1};
                relativeID(1) = relativeID(1) + 2;
            case 'timeBefore'
                timeBefore = varargin{relativeID(1)+1};
                relativeID(1) = relativeID(1) + 2;
            case 'timeAfter'
                timeAfter = varargin{relativeID(1)+1};
                relativeID(1) = relativeID(1) + 2;
            case {'Save', 'save'}
                saveOut = true;
                relativeID(1) = relativeID(1) + 1;
            case {'SaveFile', 'saveFile'}
                saveFile = varargin{relativeID(1)+1};
                relativeID(1) = relativeID(1) + 2;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{relativeID(1)});
                relativeID(1) = relativeID(1) + 1;
        end
    catch
        warning('Argument %d not recognized',relativeID(1));
        relativeID(1) = relativeID(1) + 1;
    end
end

if ~exist('ImageFiles', 'var') || isempty(ImageFiles)
    [ImageFiles,p] = uigetfile({'*.sbx;*.tif;*.imgs'}, 'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ImageFiles) % no file selected
        return
    end
    ImageFiles = fullfile(p, ImageFiles);
elseif ischar(ImageFiles)
    if isdir(ImageFiles) % directory input
        p = ImageFiles;
        ImageFiles = dir(p);
        ImageFiles = fullfile(p, {ImageFiles(~cellfun(@isempty, regexpi({ImageFiles.name}, '.*(sbx|tif)'))).name});
    else % single file input
        ImageFiles = {ImageFiles};
    end
end

if ~exist('ExperimentFile','var') || isempty(ExperimentFile)
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

if ~exist('TrialIndex', 'var') || isempty(TrialIndex)
    TrialIndex = [1 inf];
elseif islogical(TrialIndex)
    TrialIndex = find(TrialIndex);
end

if ~exist('MotionCorrect', 'var') || isempty(MotionCorrect)
    vars = whos(matfile(ExperimentFile));
    if any(strcmp(vars, 'MCdata'))
        MotionCorrect = true;
    else
        MotionCorrect = false;
    end
end


%% Determine file to save to
if saveOut && isempty(saveFile)
    if exist('ExperimentFile', 'var')
        saveFile = ExperimentFile;
    else
        saveOut = false;
    end
end


%% Load in data
if ~exist('AnalysisInfo', 'var')
    load(ExperimentFile, 'AnalysisInfo', '-mat');
end

% Determine trials to analyze
if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(1:end-1)+1:size(AnalysisInfo, 1)];
end

% Determine stimuli present in experiment
StimIDs = unique(AnalysisInfo.StimID(TrialIndex));
numStims = numel(StimIDs);


%% Motion correction
if ischar(MotionCorrect)
    load(MotionCorrect, 'MCdata', '-mat');
    MotionCorrect = true;
elseif isstruct(MotionCorrect)
    MCdata = MotionCorrect;
    MotionCorrect = true;
elseif MotionCorrect == true && ~exist('MCdata', 'var')
    load(ExperimentFile, 'MCdata', '-mat');
end
if ~exist('MCdata','var')
    MotionCorrect = false;
end

%% Load image config and determine number of frames before and after
Config = load2PConfig(ImageFiles);

numFramesBefore = round(timeBefore*Config.FrameRate/Config.Depth);
numFramesAfter = round(timeAfter*Config.FrameRate/Config.Depth);

depthID = idDepth(Config.Depth, Config.Frames, 'Depth', Depth);

%% Compute average trial
fprintf('Calculating average trials for:\n');
fprintf('\t%s\n', ImageFiles{:});

% Save settings to file
StimResponse.info.numFramesBefore = numFramesBefore;
StimResponse.info.numFramesAfter = numFramesAfter;

% Initialize outputs
numFrames = numFramesBefore+numFramesAfter+1;
AvgTrial = repmat({zeros([Config(1).size(1:2), numFrames])},numStims,1);
AvgTrialdFoF = repmat({zeros([Config(1).size(1:2), numFrames])},numStims,1);
trialdFoF = repmat({zeros([Config(1).size(1:2), numFrames])},numStims,1);
% AvgTrial = cell(numStims,1);
% AvgTrialdFoF = cell(numStims,1);
% trialdFoF = cell(numStims,1);

% Cycle through stimuli computing averages
for sindex = 1:numStims
    
    % Determine trials to average
    currentTrials = find(AnalysisInfo.StimID == StimIDs(sindex));
    currentTrials(~ismember(currentTrials, TrialIndex)) = []; % remove unwanted trials
    numTrials = numel(currentTrials);
    
    % Initialize outputs
%     numFrames = numFramesBefore+numFramesAfter+1;
%     AvgTrial{sindex} = zeros([Config(1).size(1:end-1), numFrames]);
%     AvgTrialdFoF{sindex} = zeros([Config(1).size(1:end-2), numFrames]);
%     trialdFoF{sindex} = zeros([Config(1).size(1:end-2), numTrials]);
    
    % Cycle through trials adding each to the average
    for tindex = 1:numTrials
        
        % Load trial
                StimFrames = AnalysisInfo.ExpStimFrames(currentTrials(tindex),:);
        relativeID = [find(depthID>=StimFrames(1),1,'first'), find(depthID<StimFrames(2),1,'last')];
        [frames, loadObj] = load2P(ImageFiles{1},...
            'Type',     'Direct',...
            'Depth',    Depth,...
            'Channel',  Channel,...
            'Frames',   depthID(relativeID(1)-numFramesBefore:relativeID(1)+numFramesAfter),...
            'Double');
        if MotionCorrect
            frames = applyMotionCorrection(frames, MCdata, loadObj);
        end
        frames = permute(frames,[1,2,5,3,4]);
        
        % Add trial to average
        AvgTrial{sindex} = AvgTrial{sindex} + frames; % if concerned about precision clipping high values during sum add: /numTrials;
        
        % Compute trial's dFoF
        baseline = median(frames(:,:,1:numFramesBefore),3);
        baseline(baseline<1) = 1; % in reality this never happens but don't want to enhance a small value
        frames = bsxfun(@rdivide, bsxfun(@minus, frames, baseline), baseline);
        AvgTrialdFoF{sindex} = AvgTrialdFoF{sindex} + frames; % if concerned about precision clipping high values during sum add: /numTrials;
        
        % Save for later calculation
        trialdFoF{sindex}(:,:,tindex) = mean(frames(:,:,numFramesBefore+1:numFramesBefore+diff(relativeID)+1), 3);
    end
    AvgTrial{sindex} = uint16(AvgTrial{sindex}/numTrials); % if not concerned about precision clipping high values after sum, otherwise comment out and amend above
    AvgTrialdFoF{sindex} = AvgTrialdFoF{sindex}/numTrials; % if not concerned about precision clipping high values after sum, otherwise comment out and amend above
    
    fprintf('\tfinished stim %d of %d (%d trials)\n',sindex,numStims,numTrials);
end
   

%% Compute average response
StimResponse.avg = nan([Config(1).size(1:2), numStims]);
StimResponse.se = nan([Config(1).size(1:2), numStims]);
StimResponse.pvalue = ones([Config(1).size(1:2), numStims]);
StimResponse.excited = false([Config(1).size(1:2), numStims]);
for sindex = 1:numStims
    StimResponse.avg(:,:,sindex) = mean(trialdFoF{sindex}, 3);
    StimResponse.se(:,:,sindex) = std(trialdFoF{sindex},[],3)/sqrt(size(trialdFoF{sindex},3));
    if sindex ~= 1
        [~, StimResponse.pvalue(:,:,sindex)] = ttest2(...
            trialdFoF{1},...
            trialdFoF{sindex},...
            'dim', 3);
        StimResponse.excited(:,:,sindex) = StimResponse.avg(:,:,sindex) > StimResponse.avg(:,:,1);
    end
end


%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'AvgTrial', 'AvgTrialdFoF', 'StimResponse', 'TrialIndex', '-mat', '-v7.3');
    else
        save(saveFile, 'AvgTrial', 'AvgTrialdFoF', 'StimResponse', 'TrialIndex', '-mat','-append');
    end
    fprintf('Saved average stimuli to: %s\n', saveFile);
end