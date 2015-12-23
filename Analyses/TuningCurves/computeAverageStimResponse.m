function [AvgTrial, AvgTrialdFoF, StimResponse] = computeAverageStimResponse(ImageFiles, ExperimentFile, TrialIndex, MotionCorrect, varargin)

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


%% Load in image config
Config = load2PConfig(ImageFiles);


%% Compute average trial
fprintf('Calculating average trials for:\n');
fprintf('\t%s\n', ImageFiles{:});

% Initialize outputs
AvgTrial = cell(numStims,1);
AvgTrialdFoF = cell(numStims,1);
trialdFoF = cell(numStims,1);

% Cycle through stimuli computing averages
for sindex = 1:numStims
    
    % Determine trials to average
    currentTrials = TrialIndex(AnalysisInfo.StimID(TrialIndex) == StimIDs(sindex));
    numTrials = numel(currentTrials);
    
    % Initialize outputs
    numFrames = mode(AnalysisInfo.nFrames(currentTrials));
    AvgTrial{sindex} = zeros([Config(1).size(1:end-1), numFrames]);
    AvgTrialdFoF{sindex} = zeros([Config(1).size(1:end-2), numFrames]);
    trialdFoF{sindex} = zeros([Config(1).size(1:end-2), numTrials]);
    
    % Cycle through trials adding each to the average
    for tindex = currentTrials'
        
        % Load trial
        [frames, loadObj] = load2P(AnalysisInfo.ImgFilename{tindex},...
            'Type',     'Direct',...
            'Frames',   AnalysisInfo.ExpFrames(tindex,1):AnalysisInfo.ExpFrames(tindex,1)+numFrames-1,...
            'Double');
        if MotionCorrect
            frames = applyMotionCorrection(frames, MCdata, loadObj);
        end
        
        % Add trial to average
        AvgTrial{sindex} = AvgTrial{sindex} + frames/numTrials;
        
        % Compute trial's dFoF
        lastFrame = AnalysisInfo.TrialStimFrames(tindex,1)-1;
        baseline = median(frames(:,:,:,1,1:lastFrame),5);
        baseline(baseline<1) = 1; % in reality this never happens
        frames = bsxfun(@rdivide, bsxfun(@minus, frames(:,:,:,1,:), baseline), baseline);
        AvgTrialdFoF{sindex} = AvgTrialdFoF{sindex} + permute(frames, [1,2,3,5,4])/numTrials;
        
        % Save for later calculation
        trialdFoF{sindex}(:,:,:,tindex) = mean(frames(:,:,:,1,AnalysisInfo.TrialStimFrames(tindex,1):AnalysisInfo.TrialStimFrames(tindex,2)), 5);
    end
    
    fprintf('\tfinished stim %d of %d (%d trials)\n',sindex,numStims,numTrials);
end
   

%% Compute average response
StimResponse.avg = nan([Config(1).size(1:end-2), numStims]);
StimResponse.se = nan([Config(1).size(1:end-2), numStims]);
StimResponse.pvalue = ones([Config(1).size(1:end-2), numStims]);
StimResponse.excited = false([Config(1).size(1:end-2), numStims]);
for sindex = 1:numStims
    StimResponse.avg(:,:,:,sindex) = mean(trialdFoF{sindex}, 4);
    StimResponse.se(:,:,:,sindex) = std(trialdFoF{sindex},[],4)/sqrt(size(trialdFoF{sindex},4));
    if sindex ~= 1
        [~, StimResponse.pvalue(:,:,:,sindex)] = ttest2(...
            trialdFoF{1},....
            trialdFoF{sindex},...
            'dim', 4);
        StimResponse.excited(:,:,:,sindex) = StimResponse.avg(:,:,:,sindex) > StimResponse.avg(:,:,:,1);
    end
end


%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'AvgTrial', 'AvgTrialdFoF', 'StimResponse', '-mat', '-v7.3');
    else
        save(saveFile, 'AvgTrial', 'AvgTrialdFoF', 'StimResponse', '-mat','-append');
    end
    fprintf('Saved average stimuli to: %s\n', saveFile);
end