function ROIdata = mergeROIFiles(ROIFiles, varargin)


saveOut = false;
saveFile = '';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Save'
                saveOut = true;
                index = index + 1;
            case 'SaveFile'
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

if ~exist('ROIFiles','var') || isempty(ROIFiles)
    [ROIFiles, p] = uigetfile({'*.rois;*.mat'}, 'Choose ROI file', directory,'MultiSelect','on');
    if isnumeric(ROIFiles)
        return
    end
    ROIFiles = fullfile(p,ROIFiles);
end
numFilesToAppend = numel(ROIFiles) - 1;
if ~iscell(ROIFiles) || numFilesToAppend == 0
    error('Need more than one file for merging.');
end


%% Load in data
ROIs = cell(numFilesToAppend, 1);
for findex = 1:numFilesToAppend
    load(ROIFiles{findex+1}, 'ROIdata', '-mat');
    ROIs{findex} = ROIdata;
end
load(ROIFiles{findex}, 'ROIdata', '-mat');


%% Verify data structures match
numROIs = numel(ROIdata.rois);
for findex = 1:numFilesToAppend
    if numel(ROIs{findex}.rois) ~= numROIs
        error('ROI files need the same number of ROIs');
    elseif isfield(ROIs{findex},'DataInfo') && (ROIs{findex}.DataInfo.numFramesAfter~=ROIdata.DataInfo.numFramesAfter || ROIs{findex}.DataInfo.numFramesBefore~=ROIdata.DataInfo.numFramesBefore) 
        error('ROIdata needs to be organized the same way');
    end
end


%% Merge data
for findex = 1:numFilesToAppend
    
    % Append trial info
    if isfield(ROIdata{findex},'DataInfo')
        ROIdata.DataInfo.numStimFrames = cat(1, ROIdata.DataInfo.numStimFrames, ROIs{findex}.DataInfo.numStimFrames);
        ROIdata.DataInfo.TrialIndex = cat(2, ROIdata.DataInfo.TrialIndex, ROIs{findex}.DataInfo.TrialIndex);
        ROIdata.DataInfo.StimID = cat(1, ROIdata.DataInfo.StimID, ROIs{findex}.DataInfo.StimID);
    end
    
    % Append individual roi data
    for rindex = 1:numROIs
        ROIdata.rois(rindex).rawdata = cat(2, ROIdata.rois(rindex).rawdata, ROIs{findex}.rois(rindex).rawdata);
        ROIdata.rois(rindex).rawneuropil = cat(2, ROIdata.rois(rindex).rawneuropil, ROIs{findex}.rois(rindex).rawneuropil);
        ROIdata.rois(rindex).data = cat(1, ROIdata.rois(rindex).data, ROIs{findex}.rois(rindex).data);
        ROIdata.rois(rindex).neuropil = cat(1, ROIdata.rois(rindex).neuropil, ROIs{findex}.rois(rindex).neuropil);
        ROIdata.rois(rindex).dFoF = cat(1, ROIdata.rois(rindex).dFoF, ROIs{findex}.rois(rindex).dFoF);
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

