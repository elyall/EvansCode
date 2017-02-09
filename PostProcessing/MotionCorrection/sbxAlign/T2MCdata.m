function MCdata = T2MCdata(fname, varargin)
% T2MCDATA	convert sbxalign output into MCdata structure.
%   MCdata = T2MCdata(FNAME) locates all files that fit [FNAME,'*.align'],
%   and creates MCdata structure for files that don't already have an
%   MCdata variable within them (does not access subfolders)
%
%   MCdata = T2MCdata(...,'save') saves output MCdata to its align file
%
%   MCdata = T2MCdata(...,'overwrite') recreates MCdata structure for
%   each file regardless of whether it already contains an MCdata variable

overwrite = false;
saveOut = false;

%% Parse input arguments
if ~exist('fname', 'var') || isempty(fname)
    [fname,p] = uigetfile({'*.align'},'Select align file:',directory);
    if isnumeric(fname)
        return
    end
    fname = fullfile(p, fname);
    fname = fname(1:end-8);
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

fprintf('Distributing sbxalign data for: %s\n',fname);


p = fileparts(fname);

%% Create MCdata variable

% Determine what align files exist
AlignFiles = dir([fname,'*.align']);
AlignFiles = fullfile(p,{AlignFiles(:).name});

% Determine matching sbx files
ImageFiles = closestFile(AlignFiles, '.sbx');

% Save MCdata struct to each file
temp = struct('T',{},'type',{},'date',{},'FullFilename',{},'Channel2AlignFrom',{},'Parameters',{});
if ~isempty(AlignFiles)
    for findex = 1:numel(AlignFiles)
        vars = whos(matfile(AlignFiles{findex}));
        if ~any(strcmp({vars(:).name}, 'MCdata')) || overwrite
            load(AlignFiles{findex}, 'T', '-mat');
            temp(findex).T = T;
            temp(findex).type = 'Translation';
            temp(findex).date = datestr(now);
            temp(findex).FullFilename = ImageFiles{findex};
            temp(findex).Channel2AlignFrom = 1;
            temp(findex).Parameters = [];
            if saveOut
                MCdata = temp(findex);
                save(AlignFiles{findex}, 'MCdata', '-append', '-mat');
                fprintf('\tMCdata saved to: %s\n',AlignFiles{findex});
            end
        else
            load(AlignFiles{findex},'MCdata','-mat');
            temp(findex) = MCdata;
        end
    end
end
MCdata = temp;
