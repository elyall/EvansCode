function [ROIFiles, ROIindex, FileIndex, DataSetIndex, DataSets] = loadDataSet(DictFile, MouseID)

if ~exist('MouseID', 'var')
    MouseID = 'all';
end

if ~iscell(MouseID)
    MouseID = {MouseID};
end

%% Load datasets
load(DictFile, 'DataSets', '-mat');


%% Match requests to data
StructIndex = [];
if numel(MouseID)==1 && strcmp(MouseID{1}, 'all')  % load in all datasets
    StructIndex = 1:numel(DataSets);
else
    for mindex = 1:numel(MouseID)
        temp = find(strcmp(MouseID(mindex),{DataSets(:).MouseID}));
        if ~isempty(temp)
            StructIndex = cat(2,StructIndex,temp);
        else
            warning('Dataset ''%s'' not found...', MouseID{mindex});
        end
    end
end
DataSets = DataSets(StructIndex);


%% Pull out dictionary entries
numToLoad = numel(DataSets);
ROIFiles = {};
ROIindex = [];
FileIndex = [];
DataSetIndex = [];
for mindex = 1:numToLoad
    FileIndex = cat(1, FileIndex, numel(ROIFiles)+DataSets(mindex).FileIndex);
    ROIFiles = cat(2, ROIFiles, DataSets(mindex).ROIFiles);
    ROIindex = cat(1, ROIindex, DataSets(mindex).ROIindex);
    DataSetIndex = cat(1, DataSetIndex, repmat(mindex, max(DataSets(mindex).FileIndex(:)), 1));
end
