function [ROIs, FileIndex, ROIindex] = distributeROIdata(ROIs, FieldName, Data, varargin)
% ROIs is a cell array of strings specifying the filenames of the ROI
% files, or it is a cell array of ROIdata objects
% FieldName is a field name of 'ROIdata.rois'
% Data is cell array of length numROIs

ROIindex = [];      % ROIindex is a list of ROI indices corresponding to which ROIs to pull data from
FileIndex = [];     % FileIndex is a list of file indices corresponding to which files the ROIs to pull from are in

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'ROIindex', 'ROIs', 'rois'}
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FileIndex'
                FileIndex = varargin{index+1};
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

% Convert data to cell array
if ismatrix(Data)
    Data = mat2cell(Data,ones(size(Data,1),1),size(Data,2));
elseif isnumeric(Data)
    error('Convert matrix to cell array of length numROIs');
end


%% Load ROIdata
if isempty(ROIs)
    [ROIs,p] = uigetfile({'*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p,ROIs);
end
if ~iscell(ROIs)
    ROIs = {ROIs};
end
numFiles = numel(ROIs);
if iscellstr(ROIs)
    for findex = 1:numFiles
        load(ROIs{findex}, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata;
    end
end


%% Determine ROIs to distribute data to
numROIs = numel(Data);
if isempty(ROIindex)
    numROIsPerFile = cellfun(@(x) numel(x.rois), ROIs);
    ROIindex = arrayfun(@(x) 1:x, numROIsPerFile, 'UniformOutput', false);
    ROIindex = cat(2,ROIindex{:})';
    FileIndex = repelem(1:numFiles,numROIsPerFile)';
end
if isempty(FileIndex)
    FileIndex = ones(numROIs,1); % assume all come from the first file
end
FileIDs = unique(FileIndex);
if numel(ROIindex)~=numROIs || numel(FileIndex)~=numROIs
    error('Number of Data cells to distribute is not equal to length of ROIindex/FileIndex');
end


%% Distribute data

for f = 1:numel(FileIDs)
    current = find(FileIndex==FileIDs(f));
    for r = 1:numel(current)
        ROIs{FileIDs(f)}.rois(ROIindex(current(r))).(FieldName) = Data{current(r)};
    end
end




    
    
    
        
