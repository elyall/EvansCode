function MCdata = T2MCdata(Tmats, varargin)
% T2MCDATA	convert sbxalign output into MCdata structure.
%   MCdata = T2MCdata(TMATS) creates an MCdata structure for TMATS, a T
%   matrix or a cell array of T matrices. TMATS can also contain strings
%   specifying .align files to load T from, or partial paths to find align
%   files based off of (i.e. [T,'*.align']).
%
%   MCdata = T2MCdata(...,'save') saves output MCdata to its align file if
%   its align file is specified
%
%   MCdata = T2MCdata(...,'overwrite') recreates MCdata structure for
%   each file regardless of whether it already contains an MCdata variable

ImageFile = {};
overwrite = false;
saveFile = '';
saveOut = false;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ImageFile'
                ImageFile = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'saveFile','SaveFile'}
                saveFile = varargin{index+1};
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

%% Determine files found on partial paths input
if ~iscell(Tmats)
    Tmats = {Tmats};
end

numFiles = numel(Tmats);
if ~numFiles
    MCdata = [];
    return
end

if isempty(saveFile)
    saveFile = repmat({''},numFiles,1);
elseif ~iscell(saveFile)
    saveFile = {saveFile};
end
if isempty(ImageFile)
    ImageFile = repmat({''},numFiles,1);
elseif ~iscell(ImageFile)
    ImageFile = {ImageFile};
end


%% Create MCdata variable

% Save MCdata struct to each file
temp = struct('T',{},'type',{},'date',{},'FullFilename',{},'Channel2AlignFrom',{},'Parameters',{});
for findex = 1:numFiles
    
    temp(findex).T = Tmats{findex};
    temp(findex).type = 'Translation';
    temp(findex).date = datestr(now);
    temp(findex).FullFilename = ImageFile{findex};
    temp(findex).Channel2AlignFrom = 1;
    temp(findex).Parameters = [];
    
end


%% Save to file
if saveOut && ~isempty(saveFile)
    if numel(saveFile)==numFiles
        for findex = 1:numFiles
            MCdata = temp(findex);
            save(saveFile{findex}, 'MCdata', '-append');
            fprintf('MCdata saved to: %s\n',saveFile{findex});
        end
    else
        MCdata = temp;
        save(saveFile{1}, 'MCdata', '-append');
        fprintf('MCdata saved to: %s\n',saveFile{1});
    end
end


MCdata = temp;
