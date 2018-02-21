function [newStimMat,combinations] = designStimMatrix(StimMat,InteractionTerms,ExtraVar)

[numObs,numStim] = size(StimMat);

if ~exist('InteractionTerms','var') || isempty(InteractionTerms)
    InteractionTerms = 0:numStim;
end
if ~exist('ExtraVar','var')
    ExtraVar = [];
end

% Determine combinations
combinations = {};
for index = 1:numel(InteractionTerms)
    if InteractionTerms(index)==0
        combinations = [combinations,{[]}];
    else
    currentCombs = combnk(1:numStim,InteractionTerms(index));
    if currentCombs(1)~=1 % combnk sometimes produces the output upside down from what you'd expect
        currentCombs = flipud(currentCombs);
    end
    combinations = [combinations,mat2cell(currentCombs,ones(size(currentCombs,1),1),InteractionTerms(index))'];
    end
end
numCombs = numel(combinations);

% Build stimulus matrix
newStimMat = zeros(numObs,numCombs-1); % doesn't include the zero term
for index = 1:numCombs
    if isempty(combinations{index})
        newStimMat(:,index) = ~any(StimMat,2);
    else
        newStimMat(:,index) = all(StimMat(:,combinations{index}),2);
    end
end

% Add run speed
if ~isempty(ExtraVar)
    try
        newStimMat = [newStimMat,ExtraVar];
    catch
        newStimMat = [newStimMat,ExtraVar'];
    end
end

newStimMat = double(newStimMat);
