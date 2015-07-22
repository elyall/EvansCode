function [AvgTrial, AvgTrialdFoF, StimResponse] = computeAverageStimResponse(ImageFile, ExperimentFile, TrialIndex, MotionCorrect, varargin)

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
            case 'SaveFile'
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

if ~exist('TrialIndex', 'var') || isempty(TrialIndex)
    TrialIndex = [1 inf];
end

if ~exist('MotionCorrect', 'var') || isempty(MotionCorrect)
    vars = whos(matfile(ExperimentFile));
    if any(strcmp(vars, 'MCdata'))
        MotionCorrect = true;
    else
        MotionCorrect = false;
    end
elseif isstruct(MotionCorrect)
    MCdata = MotionCorrect;
    MotionCorrect = true;
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
StimIDs = unique(AnalysisInfo.StimID);
numStims = numel(StimIDs);


%% Motion correction
if MotionCorrect == true && ~exist('MCdata', 'var')
    load(ExperimentFile, 'MCdata', '-mat');
end


%% Load in image config
Config = load2PConfig(ImageFile);


%% Compute average trial
fprintf('\nCalculating average trials for: %s\n', ImageFile);

% Initialize outputs
AvgTrial = cell(numStims,1);
AvgTrialdFoF = cell(numStims,1);
% trialdFoF = cell(numStims,1);

% Cycle through stimuli computing averages
for sindex = 1:numStims
    
    % Determine trials to average
    currentTrials = AnalysisInfo.StimID(TrialIndex) == StimIDs(sindex);
    numTrials = sum(currentTrials);
    
    % Initialize outputs
    numFrames = mode(AnalysisInfo.nFrames(currentTrials));
    AvgTrial{sindex} = zeros([Config.size(1:end-1), numFrames]);
    AvgTrialdFoF{sindex} = zeros([Config.size(1:end-2), numFrames]);
    trialdFoF{sindex} = zeros([Config.size(1:end-2), numTrials]);
    
    % Cycle through trials adding each to the average
    for tindex = find(currentTrials)'
        
        % Load trial
        [frames, loadObj] = load2P(ImageFile,...
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
StimResponse.avg = nan([Config.size(1:end-2), numStims]);
StimResponse.se = nan([Config.size(1:end-2), numStims]);
StimResponse.pvalue = ones([Config.size(1:end-2), numStims]);
StimResponse.excited = nan([Config.size(1:end-2), numStims]);
for sindex = 1:numStims
    StimResponse.avg(:,:,:,sindex) = mean(trialdFoF{sindex}, 4);
    StimResponse.se(:,:,:,sindex) = std(trialdFoF{sindex},[],4)/sqrt(size(trialdFoF{sindex},4));
    if sindex ~= 1
        [~, StimResponse.pvalue(:,:,:,sindex)] = ttest2(...
            trialdFoF{1},....
            trialdFoF{sindex},...
            'dim', 4);
        StimResponse.excited(:,:,:,sindex) = StimResponse.avg(:,:,:,sindex) >= StimResponse.avg(:,:,:,1);
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