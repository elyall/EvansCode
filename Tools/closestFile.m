function MatchedFile = closestFile(StringIn, FileExt, Dir)
%CLOSESTFILE    Determine file(s) that most closely matches input file(s).
%   MATCHEDFILE = fullfile(STRINGIN) finds the closest matching .sbx file
%   within STRINGIN's directory, by counting the number of contiguous,
%   shared characters.
%
%   MATCHEDFILE = fullfile(STRINGIN, FILEEXT) finds the closest matching
%   file ending in FILEEXT within STRINGIN's directory.
%
%   MATCHEDFILE = fullfile(STRINGIN, FILEEXT, DIR) finds the closest
%   matching file ending in FILEEXT within directory DIR. DIR can be a
%   single string or a cell array of strings of length STRINGIN


%% Parse input arguments
if ~iscell(StringIn)
    StringIn = {StringIn};
end
numFiles = numel(StringIn);

if ~exist('FileExt','var') || isempty(FileExt)
    FileExt = '.sbx';
end

if ~exist('Dir','var') || isempty(Dir)
    Dir = cellfun(@fileparts,StringIn,'UniformOutput',false);
elseif ischar(Dir)
    Dir = {Dir};
end
if numel(Dir) ~= numFiles
    Dir = repmat(Dir,numFiles,1);
end


%% Finding matching file
MatchedFile = cell(numFiles,1); % initialize output
for findex = 1:numFiles
    
    % Determine files that exist in specified directory
    Files = dir(fullfile(Dir{findex},['*',FileExt]));
    Files = {Files(:).name};
    
    % Determine closest match
    if ~isempty(Files)
        [~,fn,~] = fileparts(StringIn{findex});
        numSame = nan(numel(Files),1);
        for mindex = 1:numel(Files)
            index = 1;
            while index<=numel(fn) && index<=numel(Files{mindex}) && Files{mindex}(index)==fn(index)
                index = index + 1;
            end
            numSame(mindex) = index-1;
        end
        index = find(numSame==max(numSame));
        if numel(index)==1
            MatchedFile{findex} = fullfile(Dir{findex},Files{index});
        else
            MatchedFile{findex} = fullfile(Dir{findex},Files(index));
        end
    end
    
end

