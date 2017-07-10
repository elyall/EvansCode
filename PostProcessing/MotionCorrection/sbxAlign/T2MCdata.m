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

overwrite = false;
saveOut = false;

%% Parse input arguments
if ~exist('Tmats', 'var') || isempty(Tmats)
    [Tmats,p] = uigetfile({'*.align'},'Select align file:',directory);
    if isnumeric(Tmats)
        return
    end
    Tmats = fullfile(p,Tmats);
end
if ~iscell(Tmats)
    Tmats = {Tmats};
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case 'overwrite'
                overwrite = ~overwrite;
                index = index + 1;
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
for findex = numel(Tmats):-1:1
    if ischar(Tmats{findex}) && ~strcmp(Tmats{findex}(end-4:end),'align')
        temp = dir([Tmats{findex},'*.align']);
        temp = fullfile(temp(1).folder,{temp(:).name});
        Tmats = [Tmats,temp];
        Tmats(findex) = [];
    end
end
numFiles = numel(Tmats);
if ~numFiles
    MCdata = [];
    return
end

%% Create MCdata variable

% Save MCdata struct to each file
AlignFiles = cell(numFiles,1);
temp = struct('T',{},'type',{},'date',{},'FullFilename',{},'Channel2AlignFrom',{},'Parameters',{});
for findex = 1:numFiles
    if ischar(Tmats{findex})
        fprintf('Converting T to MCdata: %s\n',Tmats{findex});
        AlignFiles{findex} = Tmats{findex};
        ImageFile = closestFile(AlignFiles{findex}, '.sbx');
        vars = whos(matfile(Tmats{findex}));
        load(Tmats{findex}, 'T', '-mat');
        Tmats{findex} = T;
    else
        fprintf('Converting T to MCdata\n');
        ImageFile = '';
    end
    if isempty(AlignFiles{findex}) || ~any(strcmp({vars(:).name},'MCdata')) || overwrite
        temp(findex).T = Tmats{findex};
        temp(findex).type = 'Translation';
        temp(findex).date = datestr(now);
        temp(findex).FullFilename = ImageFile;
        temp(findex).Channel2AlignFrom = 1;
        temp(findex).Parameters = [];
        if saveOut && ~isempty(AlignFiles{findex})
            MCdata = temp(findex);
            save(AlignFiles{findex}, 'MCdata', '-append', '-mat');
            fprintf('\tMCdata saved to: %s\n',AlignFiles{findex});
        end
    else
        load(Tmats{findex},'MCdata','-mat');
        temp(findex) = MCdata;
    end
end
MCdata = temp;
