function [Data, FileIndex, ROIindex, Labels, Tag, ROIs] = gatherROIdata(ROIs, FieldName, varargin)
% ROIs is a cell array of strings specifying the filenames of the ROI
% files, or it is a cell array of ROIdata objects
% FieldName is a field name of 'ROIdata.rois'


ROIindex = [];      % ROIindex is a list of ROI indices corresponding to which ROIs to pull data from
FileIndex = [];     % FileIndex is a list of file indices corresponding to which files the ROIs to pull from are in
Label = false;      % string or cell array of strings of ROI labels to gather data from, or false if gathering data from ROIs with any label (use 'none' for unlabeled)
outputFormat = [];  % 'numeric' or 'cell' specifying what format the output should be in
Index = [];         % vector of linear indices specifying which data to keep
Shape = [];         % vector specifying how to shape the data for each ROI prior to concatenation

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
            case 'Index'
                Index = varargin{index+1};
                index = index + 2;
            case 'Shape'
                Shape = varargin{index+1};
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
    if ~isfield(ROIs{findex}.rois(1),'label')
        ROIs{findex}.rois(1).label = '';
    end
    unlabeledIndex = cellfun(@isempty,{ROIs{findex}.rois(:).label});
    if any(unlabeledIndex)
        [ROIs{findex}.rois(unlabeledIndex).label] = deal({'none'});
    end
end


%% Determine ROIs to gather data from
if isempty(ROIindex)
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
numFiles = numel(FileIDs);


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
for findex = 1:numFiles
    [Data{FileIndex==FileIDs(findex)}] = deal(ROIs{FileIDs(findex)}.rois(ROIindex(FileIndex==FileIDs(findex))).(FieldName)); % deal all data from current file to the cell array
end

% If cells within cells, pull out and concatenate cells
while all(cellfun(@iscell,Data)) && all(cellfun(@(x) isequal(size(Data{1}),size(x)), Data(2:end)))
    if all(cellfun(@isrow,Data))
        Data = cat(1,Data{:}); 
    elseif all(cellfun(@iscolumn,Data))
        Data = cat(2,Data{:})';
    else
        break
    end
end

% Keep only requested data & shape as desired
if ~isempty(Index)
    try
        Data = cellfun(@(x) x(Index), Data, 'UniformOutput',false);
    catch
        warning('Failed to grab requested Indices; likely Index exceeds matrix dimensions.');
    end
end
if ~isempty(Shape)
    try
        Data = cellfun(@(x) reshape(x,Shape), Data, 'UniformOutput',false);
    catch
        warning('Failed to reshape as desired; likely the number of elements does not equal the product of the dimensions');
    end
end

        
%% Convert to numeric array if desired/possible

% Convert output to matrix if possible & desired or not specified
if (strcmp(outputFormat,'numeric') || isempty(outputFormat)) && all(cellfun(@isnumeric, Data(:))) &&  all(cellfun(@(x) isequal(size(Data{1}),size(x)), Data(2:end))) % elements are all numeric and have the same size
    N = ndims(Data{1})+1;      % determine dimension along which to concatenate
    Data = cat(N,Data{:});     % concatenate data
    Data = shiftdim(Data,N-1); % move ROI dimension to the first dimension
    Data = squeeze(Data);      % remove excess dimensions
elseif strcmp(outputFormat,'numeric')
    warning('Outputting as matrix is not possible due to difference in matrix sizes or data types. Will output as cell.');
end

