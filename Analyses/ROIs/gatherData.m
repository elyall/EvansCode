function [Data, FileIndex, Index1, In] = gatherData(In, Var, Pointers, FileIndex, Index1, Index2, varargin)

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
    [In, p] = uigetfile({'*.mat'},'Choose mat file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(In)
        return
    end
    In = fullfile(p, In);
end
if ~iscell(In)
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

if ~exist('FileIndex', 'var')
    FileIndex = [];
end

if ~exist('Index1', 'var') || isempty(Index1)
    Index1 = ':';
end

if ~exist('Index2', 'var') || isempty(Index2)
    Index2 = ':';
end


%% Load data
if iscellstr(In)
    if isempty(Var)
        error('Requires specific variable name to load from files input.');
    end
    Filenames = In;
    In = cell(numFiles, 1);
    for findex = 1:numFiles
        temp = load(Filenames{findex}, Var, '-mat');
        In{findex} = temp.(Var);
    end
end

%% Pull out level desired
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


%% Determine data to pull out
dim = cellfun(@size, In, 'UniformOutput', false);
if isequal(Index1, [1,inf]) || ischar(Index1)
    numPerFile = cellfun(@(x) x(1), dim);
    Index1 = [];
    FileIndex = [];
    for findex = 1:numFiles
        Index1 = cat(2, Index1, 1:numPerFile(findex));
        FileIndex = cat(2, FileIndex, findex*ones(1,numPerFile(findex)));
    end
else
    numPerFile = zeros(numFiles, 1);
    for findex = 1:numFiles
        numPerFile(findex) = sum(FileIndex==findex);
    end
end
numTotal = sum(numPerFile);


%% Determine data size & initialize outputs
if isempty(outputFormat)
    if isnumeric(In{1}) || iscell(In{1})
        if numel(unique(cellfun(@numel,dim)))~=1 % dimensions vary across files
            outputFormat = 'cell';
        else
            outputFormat = 'numeric';
        end
    else
        outputFormat = 'cell';
    end
end


%% Pull out data
switch outputFormat
    case 'numeric'
        if ischar(Index2)
            num = size(In{1},2);
        else
            num = numel(Index2);
        end
        Data = nan(numTotal,num);
    case 'cell'
        Data = cell(numTotal,1);
end
for index = 1:numTotal
    switch outputFormat
        case 'numeric'
            Data(index,:) = In{FileIndex(index)}(Index1(index),Index2);
        case 'cell'
            Data{index} = In{FileIndex(index)}(Index1(index),Index2);
    end
end


%% Reshape data


