function [ROIindex, FileIndex, ROIs] = fixLabels(ROIFiles)

[ROIindex, FileIndex, ~, ROIs] = determineROIs(ROIFiles, 'struct-index');

for index = 1:size(ROIindex, 1)
    if isempty(ROIs{FileIndex(index, 1)}.rois(ROIindex(index, 1)).label)
        ROIs{FileIndex(index, 1)}.rois(ROIindex(index, 1)).label = {'none'};
    end
    if isempty(ROIs{FileIndex(index, 2)}.rois(ROIindex(index, 2)).label)
        ROIs{FileIndex(index, 2)}.rois(ROIindex(index, 2)).label = {'none'};
    end
    fprintf('%s \t%s\n', ROIs{FileIndex(index, 1)}.rois(ROIindex(index, 1)).label{1}, ROIs{FileIndex(index, 2)}.rois(ROIindex(index, 2)).label{1});
    if ~strcmp(ROIs{FileIndex(index, 1)}.rois(ROIindex(index, 1)).label{1}, ROIs{FileIndex(index, 2)}.rois(ROIindex(index, 2)).label{1})
        fprintf('error with ROI %d in File %d\n', ROIindex(index, 1), FileIndex(index, 1));
        ROIs{FileIndex(index, 1)}.rois(ROIindex(index, 1)).label = {'PV'};
        ROIs{FileIndex(index, 2)}.rois(ROIindex(index, 2)).label = {'PV'};
    end
end