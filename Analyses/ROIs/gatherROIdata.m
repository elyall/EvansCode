function [Data, FileIndex, ROIindex, Labels, Tag, ROIs] = gatherROIdata(ROIs, FieldName, varargin)
% ROIs is a cell array of strings specifying the filenames of the ROI
% files, or it is a cell array of ROIdata objects
% FieldName is a field name of 'ROIdata.rois'


ROIindex = 'all';   % ROIindex is a list of ROI indices corresponding to which ROIs to pull data from
FileIndex = [];     % FileIndex is a list of file indices corresponding to which files the ROIs to pull from are in
Label = false;      % string or cell array of strings of ROI labels to gather data from, or false if gathering data from ROIs with any label (use 'none' for unlabeled)
outputFormat = [];  % 'numeric' or 'cell' specifying what format the output should be in

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
            case 'Label'
                Label = varargin{index+1};
                index = index + 2;
            case 'outputFormat'
                outputFormat = varargin{index+1};
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

if ~exist('ROIs', 'var') || isempty(ROIs)
    [ROIs,p] = uigetfile({'*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p,ROIs);
end
if ~iscell(ROIs)
    ROIs = {ROIs};
end

if ~exist('FieldName', 'var') || isempty(FieldName)
    FieldName = 'dFoF'; % defaults to dFoF field
end


%% Load data
numFiles = numel(ROIs);
if iscellstr(ROIs)
    for findex = 1:numFiles
        load(ROIs{findex}, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata;
    end
end


%% Label unlabeled ROIs
for findex = 1:numFiles
    unlabeledIndex = cellfun(@isempty,{ROIs{findex}.rois(:).label});
    if any(unlabeledIndex)
        [ROIs{findex}.rois(unlabeledIndex).label] = deal({'none'});
    end
end


%% Determine ROIs to gather data from
if (isnumeric(ROIindex) && isequal(ROIindex, [1,inf])) || (ischar(ROIindex) && strcmp(ROIindex, 'all'))
    numROIsPerFile = cellfun(@(x) numel(x.rois), ROIs);
    ROIindex = arrayfun(@(x) 1:x, numROIsPerFile, 'UniformOutput', false);
    ROIindex = cat(2,ROIindex{:})';
    FileIndex = repelem(1:numFiles,numROIsPerFile)';
end
numROIs = numel(ROIindex);

if isempty(FileIndex)
    FileIndex = ones(numROIs,1); % assume all come from the first file
end
FileIDs = unique(FileIndex);
numROIsPerFile = arrayfun(@(x) nnz(FileIndex==x), FileIDs)';


%% Determine data size & initialize outputs
if isempty(outputFormat)
    temp = ROIs{FileIndex(1)}.rois(ROIindex(1)).(FieldName);
    if ~iscell(temp)
        if isequal(size(ROIs{FileIndex(1)}.rois(ROIindex(1)).(FieldName)),size(ROIs{FileIndex(end)}.rois(ROIindex(end)).(FieldName)));
            outputFormat = 'numeric';
        else
            outputFormat = 'cell';
        end
    else
        if ischar(temp) || (iscell(temp) && numel(temp)>1)
            outputFormat = 'cell';
        elseif iscell(temp)
            outputFormat = 'numeric'; % ensures not a cell within a 1x1 cell
        end
    end
end


%% Pull out data

% Pull out labels
if nargout>3 || ~isequal(Label,false)
    Labels = cell(numROIs, 1);
    Tag = cell(numROIs, 1);
    for rindex = 1:numROIs
        Labels{rindex} = ROIs{FileIndex(rindex)}.rois(ROIindex(rindex)).label{1};
        Tag{rindex} = ROIs{FileIndex(rindex)}.rois(ROIindex(rindex)).tag;
    end
end

% Remove unwanted ROIs
if ~isequal(Label,false)
    unwanted = ~ismember(Labels, Label);
    ROIindex(unwanted) = [];
    FileIndex(unwanted) = [];
    Labels(unwanted) = [];
    Tag(unwanted) = [];
end

% Grab relevant data into cell array
Data = cell(numROIs,1);
ind = [0,cumsum(numROIsPerFile)];
for findex = 1:numFiles
    [Data{ind(findex)+1:ind(findex+1)}] = deal(ROIs{FileIDs(findex)}.rois(ROIindex(FileIndex==FileIDs(findex))).(FieldName)); % deal all data from current file to the cell array
end
if strcmp(outputFormat,'numeric')
    try
        N = ndims(Data{1})+1;      % determine dimension along which to concatenate
        Data = cat(N,Data{:});     % concatenate data
        Data = shiftdim(Data,N-1); % move ROI dimension to the first dimension
    catch
        warning('Failed to output as matrix. Will output as cell.');
    end
end

