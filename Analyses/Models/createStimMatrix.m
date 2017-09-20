function StimMat = createStimMatrix(StimIndex, varargin)
% STIMINDEX is a vector of indices of length numTrials or length numFrames.
% STIMINDEX indexes into the rows of DICT, a logical matrix of size
% numConditions by numStimuli.
%
% STIMMAT is a logical matrix of size numTrials, or numFrames, by
% numStimuli. numStimuli is determined by the width of DICT, or the maximum
% index in STIMINDEX

ControlID = 0;
Dict = []; % logical numCond x numStim matrix specifying which stimuli were present during each condition

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ControlID'
                ControlID = varargin{index+1};
                index = index + 2;
            case 'Dict'
                Dict = varargin{index+1};
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


%% Set dictionary

% Format StimIndex to treat control trials the same as non-stimuli frames
StimIndex(StimIndex==ControlID) = nan; % set control trials to nan, meaning no stimulus is presented
if any(StimIndex<1)
    StimIndex = StimIndex - min(StimIndex) + 1; % ensure index is positive (assumes offsetting will fix it)
end
StimIndex(isnan(StimIndex)) = 0; % set nan frames to 0 index

% Create dictionary & ensure "no stimulus" condition is first row
if isempty(Dict)
    numStim = max(StimIndex);               % assumes stimuli exist through largest index
    Dict = [zeros(1,numStim);eye(numStim)]; % assumes each index (condition) refers to a unique stimulus
elseif any(Dict(1,:))                       % "no stimulus" condition does not exist in dictionary of conditions
    Dict = [zeros(1,size(Dict,2));Dict];    % add "no stimulus" condition as first entry (row)
end
StimIndex = StimIndex + 1; % shift index to index within Dict (i.e. account for "no stimulus" condition)


%% Create stim matrix

StimMat = Dict(StimIndex,:); % create stim matrix


