function Labels = genLabels(Labels,StimLog,Index,Label)
% Labels is name for each column in StimLog
% StimLog is numConds by numStimuli
% Index is a vector of indices into numConds
% Label is a cell array of strings that specifies the specific label for
% conditions in Index

N = size(StimLog,1);
StimLog = logical(StimLog);

Labels = arrayfun(@(x) Labels(StimLog(x,:)), 1:N, 'UniformOutput',false); % gather relevant pieces
Labels = arrayfun(@(x) sprintf('%s,',Labels{x}{:}), 1:N, 'UniformOutput',false); % concatenate text
Labels = cellfun(@(x) x(1:end-1), Labels, 'UniformOutput',false); % remove last comma

if exist('Index','var') && ~isempty(Index)
    for ind = 1:numel(Index)
        Labels{Index(ind)} = Label{ind};
    end
end