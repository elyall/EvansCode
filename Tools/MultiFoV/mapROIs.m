function [ROIs, H, W] = mapROIs(ROIFiles, ExperimentFiles, SaveFile)


if ~exist('ROIFiles', 'var')
    [ROIFiles,p] = uigetfile({'*.mat'},'Choose ROI files','MultiSelect','on');
    if isnumeric(ROIFiles)
        return
    elseif ischar(ROIFiles)
        ROIFiles = {fullfile(p,ROIFiles)};
    else
        for index = 1:numel(ROIFiles)
            ROIFiles{index} = fullfile(p, ROIFiles{index});
        end
    end
else
    [p,~,~] = fileparts(ROIFiles{1});
end

if ~exist('ExperimentFiles', 'var')
    [ExperimentFiles, p] = uigetfile({'*.mat;*.ciexp'},'Choose Experiment file',p,'MultiSelect','on');
    if isnumeric(ExperimentFiles)
        return
    end
    if isnumeric(ExperimentFiles)
        return
    elseif ischar(ExperimentFiles)
        ExperimentFiles = {fullfile(p,ExperimentFiles)};
    else
        for index = 1:numel(ExperimentFiles)
            ExperimentFiles{index} = fullfile(p, ExperimentFiles{index});
        end
    end
end

if ~exist('SaveFile', 'var')
    SaveFile = [];
end

%% Load in maps
N = numel(ExperimentFiles);
for findex = 1:N
    Data(findex) = load(ExperimentFiles{findex}, 'Map', '-mat');
end
for findex = 1:N
    Data(findex).file = ExperimentFiles{findex};
end

%% Build Map
if isempty(SaveFile)
    fprintf('Building composite map...');
    [Data, Map] = buildCompositeMap(Data);
    fprintf('\tcomplete');
else
    load(SaveFile, 'Map');
    if ~exist('Map', 'var')
        fprintf('Building composite map...');
        [Data, Map] = buildCompositeMap(Data);
        save(SaveFile, 'ExperimentFiles', 'Map', '-mat', '-append');
        fprintf('\tcomplete');
    else
        Data = buildCompositeMap(Data);
    end
end
[H, W, ~] = size(Map);

%% Load in ROIs
nrois = 0;
for findex = 1:N
    dy = Data(findex).ylim(1);
    dx = Data(findex).xlim(1);
    load(ROIFiles{findex},'ROIdata', '-mat');
    for rindex = 1:numel(ROIdata.rois)
        ROIs.rois(nrois+rindex) = ROIdata.rois(rindex);
        ROIs.rois(nrois+rindex).vertices(:,2) = ROIs.rois(nrois+rindex).vertices(:,2) + dy;
        ROIs.rois(nrois+rindex).vertices(:,1) = ROIs.rois(nrois+rindex).vertices(:,1) + dx;
    end
    nrois = nrois + numel(ROIdata.rois);
end

