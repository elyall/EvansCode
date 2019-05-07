function [AvgTrial, AvgTrialdFoF, StimResponse] = computeAverageStimResponse2(ImageFiles, AnalysisInfo, MCdata, varargin)

Depth = 1;
Channel = 1;
timeBefore = 1;
timeAfter = 3;
% Filter = [];
Filter = fspecial('gaussian',5,1);
Baseline = 'catch'; % scalar specifying the trial-wise prctile, 'trial-wise', 'catch', or a frame that represents the baseline
TrialIndex = [1,inf];
numFramesBaseline = 16;

saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'TrialIndex','Trials','trials'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case {'Depth','depth'}
                Depth = varargin{index+1};
                index = index + 2;
            case 'timeBefore'
                timeBefore = varargin{index+1};
                index = index + 2;
            case 'timeAfter'
                timeAfter = varargin{index+1};
                index = index + 2;
            case {'filter','Filter'}
                Filter = varargin{index+1};
                index = index + 2;
            case 'Baseline'
                Baseline = varargin{index+1};
                index = index + 2;
            case 'numFramesBaseline'
                numFramesBaseline = varargin{index+1};
                index = index + 2;
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

if ~exist('ImageFiles', 'var') || isempty(ImageFiles)
    [ImageFiles,p] = uigetfile({'*.sbx;*.tif;*.imgs'}, 'Choose files to process', directory, 'MultiSelect', 'on');
    if isnumeric(ImageFiles) % no file selected
        return
    end
    ImageFiles = fullfile(p, ImageFiles);
elseif ischar(ImageFiles) && isdir(ImageFiles) % directory input
        p = ImageFiles;
        ImageFiles = dir(p);
        ImageFiles = fullfile(p, {ImageFiles(~cellfun(@isempty, regexpi({ImageFiles.name}, '.*(sbx|tif)'))).name});
end
if ischar(ImageFiles) % single file input
    ImageFiles = {ImageFiles};
end


%% Load in data

% Experiment info
if ~exist('AnalysisInfo','var') || isempty(AnalysisInfo)
    AnalysisInfo = closestFile(ImageFiles{1},'.exp'); % guess experiment file
    AnalysisInfo = AnalysisInfo{1};
    if isempty(AnalysisInfo)
        [AnalysisInfo, p] = uigetfile({'*.mat'},'Choose Experiment file',directory);
        if isnumeric(AnalysisInfo)
            return
        end
        AnalysisInfo = fullfile(p,AnalysisInfo);
    end
end
if ischar(AnalysisInfo)
    load(AnalysisInfo, 'AnalysisInfo', '-mat');
end

% Motion correction data
if ~exist('MCdata','var') || isempty(MCdata)
    MotionCorrect = false;
else
    MotionCorrect = true;
    if ischar(MCdata)
        load(MCdata, 'MCdata', '-mat');
    end
end


%% Determine analysis parameters

% Determine trials to analyze
if islogical(TrialIndex)
    TrialIndex = find(TrialIndex);
end
if TrialIndex(end) == inf
    TrialIndex = [TrialIndex(1:end-1), TrialIndex(1:end-1)+1:size(AnalysisInfo, 1)];
end

% Determine stimuli present in experiment
StimIDs = unique(AnalysisInfo.StimID(TrialIndex));
numStims = numel(StimIDs);

% Determine number of frames to gather and which frames correspond to the
% current depth
Config = load2PConfig(ImageFiles);
numFramesBefore = round(timeBefore*Config.FrameRate/Config.Depth);
numFramesAfter = round(timeAfter*Config.FrameRate/Config.Depth);
depthID = idDepth(Config.Depth, Config.Frames, 'Depth', Depth);


%% Compute baseline

switch Baseline
    case 'catch'
        Index = AnalysisInfo.StimID==0; % index of catch stimuli
        Index = logical([0;Index(1:end-1)]); % index of stimuli after catch
        frames = AnalysisInfo.ExpStimFrames(Index,1)-1; % last frame before following stimulus
        frames = [max(frames-numFramesBaseline+1,1),frames]; % list of first and last frames to average
        FrameIndex = false(Config.Frames,1); % index of frames to average
        for t = 1:size(frames,1)
            FrameIndex(frames(t,1):frames(t,2)) = true;
        end
        currentDepth = false(Config.Frames,1);
        currentDepth(Depth:Config.Depth:end) = true; % index of frames corresponding to current depth
        FrameIndex = all([FrameIndex,currentDepth],2);
        
        [frames, ~] = load2P(ImageFiles,...
            'Type',     'Direct',...
            'Depth',    Depth,...
            'Channel',  Channel,...
            'Frames',   find(FrameIndex),...
            'Double');
        Baseline = mean(frames,5);
end


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
totalFrames = Config.Frames;
ExpStimFrames = AnalysisInfo.ExpStimFrames;
StimID = AnalysisInfo.StimID;
for sindex = 1:numStims
    
    % Determine trials to average
    currentTrials = find(StimID == StimIDs(sindex));
    currentTrials(~ismember(currentTrials, TrialIndex)) = []; % remove unwanted trials
    numTrials = numel(currentTrials);
    
    % Initialize outputs
%     numFrames = numFramesBefore+numFramesAfter+1;
%     AvgTrial{sindex} = zeros([Config(1).size(1:end-1), numFrames]);
%     AvgTrialdFoF{sindex} = zeros([Config(1).size(1:end-2), numFrames]);
%     trialdFoF{sindex} = zeros([Config(1).size(1:end-2), numTrials]);
    
    % Cycle through trials adding each to the average
    badTrials = 0;
    for tindex = 1:numTrials
        
        % Load trial
        StimFrames = ExpStimFrames(currentTrials(tindex),:);
        relativeID = [find(depthID>=StimFrames(1),1,'first'), find(depthID<StimFrames(2),1,'last')];
        if relativeID(1)+numFramesAfter > totalFrames % trial is last trial of experiment and has less frames available than requested
            badTrials = badTrials+1; % should occur at most once per experiment (and only if experiment was stopped early or timeAfter is large)
            continue
        end
        [frames, loadObj] = load2P(ImageFiles,...
            'Type',     'Direct',...
            'Depth',    Depth,...
            'Channel',  Channel,...
            'Frames',   depthID(relativeID(1)-numFramesBefore:relativeID(1)+numFramesAfter),...
            'Double');
        if ~isempty(Filter)
            frames = imfilter(frames, Filter);
        end
        if MotionCorrect
            frames = applyMotionCorrection(frames, MCdata, loadObj);
        end
        frames = permute(frames,[1,2,5,3,4]);
        
        % Add trial to average
        AvgTrial{sindex} = AvgTrial{sindex} + frames; % if concerned about precision clipping high values during sum add: /numTrials;
        
        % Compute trial's dFoF
        if ischar(Baseline) && strcmp(Baseline,'trial-wise')
            baseline = median(frames(:,:,1:numFramesBefore),3);
        elseif isscalar(Baseline)
            baseline = prctile(frames,Baseline,3);
        else
            baseline = Baseline;
        end
        baseline(baseline<1) = 1; % in reality this never happens but don't want to enhance a small value
        frames = bsxfun(@rdivide, bsxfun(@minus, frames, baseline), baseline);
        AvgTrialdFoF{sindex} = AvgTrialdFoF{sindex} + frames; % if concerned about precision clipping high values during sum, add: /numTrials;
        
        % Save for later calculation
        trialdFoF{sindex}(:,:,tindex) = mean(frames(:,:,numFramesBefore+1:numFramesBefore+diff(relativeID)+1), 3);
    end
    numTrials = numTrials - badTrials;
    AvgTrial{sindex} = uint16(AvgTrial{sindex}/numTrials); % if concerned about precision clipping high values after sum, comment out and amend above
    AvgTrialdFoF{sindex} = AvgTrialdFoF{sindex}/numTrials; % if concerned about precision clipping high values after sum, comment out and amend above
    
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
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile, 'file')
        save(saveFile, 'AvgTrial', 'AvgTrialdFoF', 'StimResponse', 'TrialIndex', '-mat', '-v7.3');
    else
        save(saveFile, 'AvgTrial', 'AvgTrialdFoF', 'StimResponse', 'TrialIndex', '-mat','-append');
    end
    fprintf('Saved average stimuli to: %s\n', saveFile);
end