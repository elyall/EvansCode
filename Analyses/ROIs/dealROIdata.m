function ROIdata = dealROIdata(ROIdata, Field, Data, varargin)

ROIindex = [1 inf];

saveOut = false;
saveFile = '';

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case {'Save', 'save'}
                saveOut = true;
                index = index + 1;
            case {'SaveFile', 'saveFile'}
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

if ~exist('ROIdata','var') || isempty(ROIdata)
    [ROIdata, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end


%% Load in ROIdata
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
end
if saveOut && isempty(saveFile)
    saveOut = false;
end

% Determine ROIs to distribute to
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIindex);

% Validate dimensions of data
if size(Data,1) ~= numROIs
    error('First dimension of data input has to have length equal to number of ROIs being dealt.');
end


%% Distribute data to ROIdata
if ~iscell(Data)
    sz = size(Data);
    Data = reshape(Data, sz(1), prod(sz(2:end)));
    for rindex = 1:numROIs
        if ~ismatrix(sz)
            ROIdata.rois(ROIindex(rindex)).(Field) = reshape(Data(1,:), sz(2:end));
        else
            ROIdata.rois(ROIindex(rindex)).(Field) = Data(rindex,:);
        end
    end
else
    for rindex = 1:numROIs
        ROIdata.rois(ROIindex(rindex)).(Field) = Data{rindex};
    end
end


%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end

