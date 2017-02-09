function ROIs = segment2ROIdata(fname, varargin)
% SEGMENT2ROIDATA	convert sbxsegment output into ROIdata structure.
%   ROIS = segment2ROIdata(FNAME) locates all files that fit
%   [FNAME,'*.segment'], and creates ROIdata structure for files that don't
%   already have a corresponding '.rois' file with ROIdata in it (does not
%   access subfolders)
%
%   ROIS = segment2ROIdata(...,'save') saves output ROIdata to a mat file
%   with the same name as the segment file but with a '.rois' extension
%
%   ROIS = segment2ROIdata(...,'overwrite') recreates ROIdata structure for
%   each file regardless of whether a corresponding ROIdata structure
%   already exists

overwrite = false;
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

fprintf('Distributing sbxsegment data for: %s\n',fname);


p = fileparts(fname);

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
        while ~isnan(str2double(SegmentFiles{findex}(index+5+temp+1)))
            temp = temp + 1;
        end
        Depth(findex) = str2double(SegmentFiles{findex}(index+5:index+5+temp));
    end
end

% Determine matching sbx files
ImageFiles = closestFile(SegmentFiles, '.sbx');

% Create ROIdata
ROIs = cell(numFiles,1);
if ~isempty(SegmentFiles)
    for findex = 1:numFiles
        ROIFile = [SegmentFiles{findex}(1:end-7),'rois'];
        vars = whos(matfile(ROIFile));
        if ~any(strcmp({vars(:).name}, 'ROIdata')) || overwrite
            
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
            
        else
            load(ROIFile,'ROIdata','-mat')
            ROIs{findex} = ROIdata;
        end
    end
end

% Maintain legacy output
if numel(ROIs)==1
    ROIs = ROIs{1};
end

