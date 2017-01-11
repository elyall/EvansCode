function [Data, FileIndex, ROIindex, In] = gatherData(In, Var, Pointers, Index, ROIindex, FileIndex, varargin)
% ROIs is a cell array of strings specifying the filenames of the ROI
% files, or it is a cell array of ROIdata objects
% FieldName is a field name of 'ROIdata.rois'
% Index is a number indexing in that field (e.g. '1' or use '':'' for all elements)
% Label is a cell array of strings of ROIs to gather data from, or false if
% gathering data from ROIs with any label (use 'none' for unlabeled)
% ROIindex is a list of roi indices, correspondings to the files listed in
% FileIndex

outputFormat = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
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

if ~exist('In', 'var') || isempty(In)
    directory = cd;
    [In, p] = uigetfile({'*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(In)
        return
    elseif iscell(In)
        for findex = 1:numel(In)
            In{findex} = fullfile(p, In{findex});
        end
    elseif ischar(In)
        In = {fullfile(p, In)};
    end
elseif ~iscell(In)
    In = {In};
end
numFiles = numel(In);

if ~exist('Var','var')
    Var = [];
end

if ~exist('Pointers', 'var')
    Pointers = {[]};
elseif ~iscell(Pointers)
    Pointers = {Pointers};
end
if numel(Pointers) < numFiles
    Pointers = repmat(Pointers,numFiles,1);
end

if ~exist('Index', 'var') || isempty(Index)
    Index = ':';
end

if ~exist('ROIindex', 'var') || isempty(ROIindex)
    ROIindex = 'all';
end

if ~exist('FileIndex', 'var')
    FileIndex = [];
end

%% Load data
if iscellstr(In)
    if isempty(Var)
        error('Requires specific variable name to load from files input.');
    end
    Files = In;
    In = cell(numFiles, 1);
    for findex = 1:numFiles
        temp = load(Files{findex}, Var, '-mat');
        In{findex} = temp.(Var);
    end
end

%% Organize data
for findex = 1:numFiles
    for index = 1:numel(Pointers{findex})
        if isa(Pointers{findex},'char')
            if numel(In{findex})~=1
                In{findex} = cat(1,In{findex}(:).(Pointers{index}));
            else
                In{findex} = In{findex}.(Pointers{index});
            end
        elseif isa(Pointers{findex},'numeric')
            if isa(In{findex},'cell')
                In{findex} = In{findex}(Pointers{index});
            elseif isa(In{findex},'cell')
                In{findex} = In{findex}{Pointers{index}};
            elseif isa(In{findex},'numeric')
                In{findex} = In{findex}(Pointers{index},:);
            end
        end
    end
end


%% Determine ROIs to gather data from
if (isnumeric(ROIindex) && isequal(ROIindex, [1,inf])) || (ischar(ROIindex) && strcmp(ROIindex, 'all'))
    numROIs = cellfun(@size, In, 'UniformOutput', false);
    numROIs = cellfun(@(x) x(1), numROIs);
    ROIindex = [];
    FileIndex = [];
    for findex = 1:numFiles
        ROIindex = cat(2, ROIindex, 1:numROIs(findex));
        FileIndex = cat(2, FileIndex, findex*ones(1,numROIs(findex)));
    end
else
    numROIs = zeros(numFiles, 1);
    for findex = 1:numFiles
        numROIs(findex) = sum(FileIndex==findex);
    end
end
totalROIs = sum(numROIs);


%% Determine data size & initialize outputs
if isempty(outputFormat)
    dim = cellfun(@size, In, 'UniformOutput', false);

    if numel(unique(cellfun(@numel,dim)))~=1
        outputFormat = 'cell';
    else
%         dim = cellfun(@(x) x(2:end), dim);
        outputFormat = 'numeric';
    end
end


%% Pull out data
switch outputFormat
    case 'numeric'
        if ischar(Index)
            num = size(In{1},2);
        else
            num = numel(Index);
        end
        Data = nan(totalROIs, num);
    case 'cell'
        Data = cell(totalROIs, 1);
end
for rindex = 1:totalROIs
    switch outputFormat
        case 'numeric'
            Data(rindex,:) = In{FileIndex(rindex)}(ROIindex(rindex),Index);
        case 'cell'
            Data{rindex} = In{FileIndex(rindex)}(ROIindex(rindex,Index));
    end
end


%% Reshape data


