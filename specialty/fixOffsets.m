function fixOffsets(ROIFiles)


if ~exist('ROIFiles','var') || isempty(ROIFiles)
    directory = CanalSettings('DataDirectory');
    [ROIFiles, p] = uigetfile({'*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(ROIFiles)
        return
    end
    if iscell(ROIFiles)
        for findex = 1:numel(ROIFiles)
            ROIFiles{findex} = fullfile(p,ROIFiles{findex});
        end
    elseif ischar(ROIFiles)
        ROIFiles = {fullfile(p,ROIFiles)};
    end
elseif ischar(ROIFiles)
    ROIFiles = {fullfile(p,ROIFiles)};
end

for findex = 1:numel(ROIFiles)
    load(ROIFiles{findex}, 'offset');
    offset = [0,0];
    load(ROIFiles{findex}, 'ROIdata', '-append');
end