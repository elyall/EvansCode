function saveFile = saveDataSet(saveFile, MouseID, varargin)


%% Load file's contents and determine if dataset already exists
if exist(saveFile, 'file')
    load(saveFile, 'DataSets', '-mat');
end
if ~exist('DataSets', 'var') % previous file does not exist or it does not contain any datasets
    dindex = 1;
    DataSets = struct('MouseID', [], 'ROIFiles', {}, 'ROIindex', [], 'FileIndex', []);
    DataSets(dindex).MouseID = MouseID;
else
    dindex = find(strcmp(MouseID,{DataSets(:).MouseID}));
    if isempty(dindex) % dataset doesn't previously exist
        dindex = numel(DataSets)+1;
        DataSets(dindex).MouseID = MouseID;
    end
end


%% Append new data
iindex = 1;
while iindex < numel(varargin)
%     if isstruct(varargin{iindex})
%         for fn = fieldnames(varargin{iindex})'
%             DataSets(dindex).(fn{1}) = varargin{iindex}.(fn{1});
%         end
%         iindex = iindex + 1;
    if ismember(varargin{iindex}, {'ROIFiles', 'ROIindex', 'FileIndex', 'UserData', 'PWCZ', 'DistBtwn', 'insidePWC', 'TrialIndex'})
        DataSets(dindex).(varargin{iindex}) = varargin{iindex+1};
        iindex = iindex + 2;
    end
end


%% Save to file
if exist(saveFile, 'file')
    save(saveFile, 'DataSets', '-mat', '-append');
else
    save(saveFile, 'DataSets', '-mat', '-v7.3');
end
fprintf('Saved data to %s in %s\n', MouseID, saveFile);

