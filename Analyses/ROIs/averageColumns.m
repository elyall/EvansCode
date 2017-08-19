function out = averageColumns(Data,Group,StimIndex,varargin)

Layout = [0,3,0;1,2,5;0,4,0];

% Determine parameters of data
[numROIs,numTrials,numFrames] = size(Data);
Dim = size(Layout);
[StimIDs,~,StimIndex] = unique(StimIndex);
numStim = numel(StimIDs);
[GroupIDs,~,Group] = unique(Group);
if any(GroupIDs==0)
    GroupIDs(1) = [];
    Group = Group - 1;
end
numGroups = numel(GroupIDs);

% Place outputs accordingly
out = zeros(Dim(1),Dim(2),numFrames,numStim);
for gindex = 1:numGroups
    [y,x] = find(Layout==GroupIDs(gindex));
    for sindex = 1:numStim
        out(y,x,:,sindex) = mean(mean(Data(Group==gindex,StimIndex==sindex,:),2),1); % mean across trials then across ROIs
    end
end
