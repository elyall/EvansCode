function [ROIs,MCdata] = sbxDistribute(fname, varargin)

override = false;

saveOut = false;

%% Parse input arguments
if ~exist('fname', 'var') || isempty(fname)
    [fname,p] = uigetfile({'*.segment'},'Select segment file:',directory);
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
            case 'override'
                override = ~override;
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

fprintf('Distributing sbxalign and sbxsegment data for: %s\n',fname);


p = fileparts(fname);

%% Create MCdata variable

% Determine what align files exist
AlignFiles = dir([fname,'*.align']);
AlignFiles = fullfile(p,{AlignFiles(:).name});

% Determine matching sbx files
ImageFiles = cell(numel(AlignFiles),1);
temp = dir(fullfile(p,'*.sbx'));
temp = {temp(:).name};
for findex = 1:numel(AlignFiles)
    [~,fn,~] = fileparts(AlignFiles{findex});
    numSame = nan(numel(temp),1);
    for mindex = 1:numel(temp)
        index = 1;
        while index<=numel(fn) && index<=numel(temp{mindex}) && temp{mindex}(index)==fn(index)
            index = index + 1;
        end
        numSame(mindex) = index-1;
    end
    [~,index] = max(numSame);
    ImageFiles{findex} = fullfile(p,temp{index});
end

% Save MCdata struct to each file
temp = struct('type',num2cell(nan(numel(AlignFiles),1)));
if ~isempty(AlignFiles)
    for findex = 1:numel(AlignFiles)
        vars = whos(matfile(AlignFiles{findex}));
        if ~any(strcmp({vars(:).name}, 'MCdata')) || override
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
        end
    end
end
MCdata = temp;


%% Create ROIdata

% Determine what segment files exist
SegmentFiles = dir([fname,'*.segment']);
SegmentFiles = fullfile(p,{SegmentFiles(:).name});
numFiles = numel(SegmentFiles);

% Determine depth of segment files
Depth = nan(numFiles,1);
for findex = 1:numFiles
    index = strfind(SegmentFiles{findex},'depth');
    if isempty(index)
        Depth(findex) = 1;
    else
        temp = 0;
        while ~isnan(str2double(SegmentFiles{findex}(index+5+1)))
            temp = temp + 1;
        end
        Depth(findex) = str2double(SegmentFiles{findex}(index+5:index+5+temp));
    end
end

% Determine matching sbx files
ImageFiles = cell(numFiles,1);
temp = dir(fullfile(p,'*.sbx'));
temp = {temp(:).name};
for findex = 1:numFiles
    [~,fn,~] = fileparts(SegmentFiles{findex});
    numSame = nan(numel(temp),1);
    for mindex = 1:numel(temp)
        index = 1;
        while index<=numel(fn) && index<=numel(temp{mindex}) && temp{mindex}(index)==fn(index)
            index = index + 1;
        end
        numSame(mindex) = index-1;
    end
    [~,index] = max(numSame);
    ImageFiles{findex} = fullfile(p,temp{index});
end

% Create ROIdata
ROIs = cell(numFiles,1);
if ~isempty(SegmentFiles)
    for findex = 1:numFiles
        ROIFile = [SegmentFiles{findex}(1:end-7),'rois'];
        vars = whos(matfile(ROIFile));
        if ~any(strcmp({vars(:).name}, 'ROIdata')) || override
            
            % Create ROIdata
            ROIs{findex} = createROIdata(SegmentFiles{findex}, 'ImageFile', ImageFiles(findex), 'Depth', Depth(findex));
            
            % Save ROIdata to file
            if saveOut
                ROIdata = ROIs{findex};
                if ~exist(ROIFile, 'file')
                    save(ROIFile, 'ROIdata', '-mat', '-v7.3');
                else
                    save(ROIFile, 'ROIdata', '-mat', '-append');
                end
                fprintf('\tROIdata saved to: %s\n',ROIFile);
            end
            
        end
    end
end

% Maintain legacy output
if numel(ROIs)==1
    ROIs = ROIs{1};
end

