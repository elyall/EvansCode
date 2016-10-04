function ROIs = computeCentroid(ROIs, ROIFiles)

saveOut = true;

%% Load data
numFiles = numel(ROIs);
if iscellstr(ROIs)
    ROIFiles = ROIs;
    ROIs = cell(numFiles, 1);
    for findex = 1:numFiles
        load(ROIFiles, 'ROIdata');
        ROIs{findex} = ROIdata;
    end
end

%% Compute centroid
for findex = 1:numFiles
    for rindex = 1:numel(ROIs{findex}.rois)
        temp = regionprops(ROIs{findex}.rois(rindex).pixels, 'centroid');
        ROIs{findex}.rois(rindex).centroid = temp.Centroid;
    end
end

%% Save to ROI file
if saveOut && exist('ROIFiles', 'var')
    for findex = 1:numFiles
        ROIdata = ROIs{findex};
        save(ROIFiles{findex}, 'ROIdata', '-append', '-mat');
        fprintf('\nROIdata saved to: %s', ROIFiles{findex});
    end
end