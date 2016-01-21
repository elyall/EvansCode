function [ROIFiles, ROIindex, FileIndex, PWCZ, DataSets] = loadDataSet(DictFile, MouseID)

if ~exist('MouseID', 'var')
    MouseID = 'all';
end

if ~iscell(MouseID)
    MouseID = {MouseID};
end

%% Load datasets
load(DictFile, 'DataSets', '-mat');


%% Match requests to data
if numel(MouseID)==1 && strcmp(MouseID{1}, 'all')  % load in all datasets
    MouseID = num2cell(1:numel(DataSets));
else
    for mindex = numel(MouseID):-1:1
        temp = find(strcmp(MouseID(mindex),{DataSets(:).MouseID}));
        if ~isempty(temp)
            MouseID{mindex} = temp;
        else
            MouseID(mindex) = [];
        end
    end
end


%% Pull out dictionary entries
numToLoad = numel(MouseID);
ROIFiles = {};
ROIindex = [];
FileIndex = [];
PWCZ = {};
for mindex = 1:numToLoad
    FileIndex = cat(1, FileIndex, numel(ROIFiles)+DataSets(MouseID{mindex}).FileIndex);
    ROIFiles = [ROIFiles, DataSets(MouseID{mindex}).ROIFiles];
    ROIindex = cat(1, ROIindex, DataSets(MouseID{mindex}).ROIindex);
    PWCZ = [PWCZ; repmat({DataSets(MouseID{mindex}).PWCZ}, max(DataSets(MouseID{mindex}).FileIndex(:)), 1)];
end

