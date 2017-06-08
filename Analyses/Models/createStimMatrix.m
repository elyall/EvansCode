function [stim,keptFrames] = createStimMatrix(StimIndex, Stimulus, Trial, TrialIndex)
%
% StimIndex - cell array of strings with stim combinations or logical
% matrix specifying which stimulus was presented during that stimulus
%
% Stimulus - vector of NaNs w/ StimIDs during frames in which stim was
% presented
%
% Trial - vector of NaNs w/ TrialIDs for frames that fell within that trial
%
% TrialIndex - vector of trial IDs to keep

if ~exist('Trial','var')
    Trial = [];
end
if ~exist('TrialIndex','var') || isempty(TrialIndex)
    TrialIndex = [1 inf];
end

%% Format stim list
if ~logical(StimIndex)
    if isnumeric(StimIndex)                                                          % matrix input
        StimIndex = mat2cell(StimIndex,ones(size(StimIndex,1),1),size(StimIndex,2)); % assume each row is a different stimulus
    elseif iscellstr(StimIndex)                                                      % cell array of strings input
        StimIndex = cellfun(@str2num,StimIndex,'UniformOutput',false);               % convert to cell array of vectors
    end
    StimList = unique([StimIndex{:}]); % determine condition ID #s
    numStim = numel(StimIndex);        % determine # of unique stimuli
    StimCond = StimIndex;
    StimIndex = false(numStim,numel(StimList));
    for index = 1:numStim
        StimIndex(index,ismember(StimList,StimCond{index})) = true; % set conditions in stimulus to true
    end
end
[numStim,numCond] = size(StimIndex);

%% Format stimulus clock
if any(Stimulus==0)
    Stimulus = Stimulus+1;
end
Stimulus(isnan(Stimulus)) = 0;
if ~isempty(Trial)
    if TrialIndex(end)==inf
        TrialIndex = cat(2,TrialIndex(1:end-1),TrialIndex(end-1)+1:max(Trial));
    end
    keptFrames = ismember(Trial,TrialIndex);
    Stimulus = Stimulus(keptFrames);
end
numSamples = numel(Stimulus);

%% Create stim matrix
stim = zeros(numCond,numSamples);
for index = 1:numStim
    stim(:,Stimulus==index) = repmat(StimIndex(index,:)',1,nnz(Stimulus==index));
end
