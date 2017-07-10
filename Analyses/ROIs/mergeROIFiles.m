function ROIdata = mergeROIFiles(ROIs, varargin)

ROIdict = []; % numROIs in first file x numFiles to append

saveOut = false;
saveFile = '';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIdict'
                ROIdict = varargin{index+1};
                index = index + 2;
            case {'save','Save'}
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

if ~exist('ROIs','var') || isempty(ROIs)
    [ROIs, p] = uigetfile({'*.rois;*.mat'}, 'Choose ROI file', directory,'MultiSelect','on');
    if isnumeric(ROIs)
        return
    end
    ROIs = fullfile(p,ROIs);
end
numFiles = numel(ROIs);
if ~iscell(ROIs) || numFiles <= 1
    error('Need more than one file for merging.');
end


%% Load in data
if iscellstr(ROIs)
    for findex = 1:numFiles
        load(ROIs{findex}, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata;
    end
end

% Determine ROIs to append
numROIs = cellfun(@(x) numel(x.rois), ROIs);
if isempty(ROIdict)
    ROIdict = nan(numROIs(1),numFiles);
    for findex = 1:numFiles-1
        ROIdict(1:min(numROIs([1,findex])),findex) = 1:min(numROIs([1,findex+1]));
    end
end


%% Verify data structures match
for findex = 2:numFiles
    if isfield(ROIs{findex},'DataInfo') && (ROIs{findex}.DataInfo.numFramesAfter~=ROIs{1}.DataInfo.numFramesAfter || ROIs{findex}.DataInfo.numFramesBefore~=ROIs{1}.DataInfo.numFramesBefore) 
        error('ROIdatas needs to be organized the same way');
    end
end


%% Merge data
for findex = 2:numFiles
    
    % Append trial info
    ROIs{1}.DataInfo.numStimFrames = cat(1, ROIs{1}.DataInfo.numStimFrames, ROIs{findex}.DataInfo.numStimFrames);
    ROIs{1}.DataInfo.TrialIndex = cat(2, ROIs{1}.DataInfo.TrialIndex, ROIs{findex}.DataInfo.TrialIndex + max(ROIs{1}.DataInfo.TrialIndex));
    ROIs{1}.DataInfo.StimID = cat(1, ROIs{1}.DataInfo.StimID, ROIs{findex}.DataInfo.StimID);
    
    % Append individual roi data
    for rindex = 1:numROIs(1)
        if ~isnan(ROIdict(rindex,findex-1))
            index = ROIdict(rindex,findex-1);
            ROIs{1}.rois(rindex).rawdata = cat(2, ROIs{1}.rois(rindex).rawdata, ROIs{findex}.rois(index).rawdata);
            ROIs{1}.rois(rindex).rawneuropil = cat(2, ROIs{1}.rois(rindex).rawneuropil, ROIs{findex}.rois(index).rawneuropil);
            ROIs{1}.rois(rindex).data = cat(1, ROIs{1}.rois(rindex).data, ROIs{findex}.rois(index).data);
            ROIs{1}.rois(rindex).neuropil = cat(1, ROIs{1}.rois(rindex).neuropil, ROIs{findex}.rois(index).neuropil);
            ROIs{1}.rois(rindex).dFoF = cat(1, ROIs{1}.rois(rindex).dFoF, ROIs{findex}.rois(index).dFoF);
            ROIs{1}.rois(rindex).stimMean = cat(1, ROIs{1}.rois(rindex).stimMean, ROIs{findex}.rois(index).stimMean);
        else % fill with nans
            ROIs{1}.rois(rindex).rawdata = cat(2, ROIs{1}.rois(rindex).rawdata, nan(size(ROIs{findex}.rois(1).rawdata)));
            ROIs{1}.rois(rindex).rawneuropil = cat(2, ROIs{1}.rois(rindex).rawneuropil, nan(size(ROIs{findex}.rois(1).rawneuropil)));
            ROIs{1}.rois(rindex).data = cat(1, ROIs{1}.rois(rindex).data, nan(size(ROIs{findex}.rois(1).data)));
            ROIs{1}.rois(rindex).neuropil = cat(1, ROIs{1}.rois(rindex).neuropil, nan(size(ROIs{findex}.rois(1).neuropil)));
            ROIs{1}.rois(rindex).dFoF = cat(1, ROIs{1}.rois(rindex).dFoF, nan(size(ROIs{findex}.rois(1).dFoF)));
            ROIs{1}.rois(rindex).stimMean = cat(1, ROIs{1}.rois(rindex).stimMean, nan(size(ROIs{findex}.rois(1).stimMean)));
        end
    end
end
ROIdata = ROIs{1};


%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end

