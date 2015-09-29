function [TMI, ROIs] = computeTMI(ROIs, saveOut, ROIFiles)
% assumes files are input in order of [pre, post, pre, post, etc]

pos1 = 1;
pos2 = 8;
control = 1; % 1 or 0

%% Load data
numFiles = numel(ROIs);
fprintf('computeTMI - %d files:', numFiles);
if iscellstr(ROIs)
    fprintf('\tLoading files...');
    ROIFiles = ROIs;
    ROIs = cell(numFiles, 1);
    for findex = 1:numFiles
        load(ROIFiles{findex}, 'ROIdata');
        ROIs{findex} = ROIdata;
    end
end
numROIs = cellfun(@(x) numel(x.rois), ROIs);


%% Compute TMI
fprintf('\tComputing TMI...');
TMI = cell(numFiles/2, 1);
for findex = 1:2:numFiles
    TMI{(findex+1)/2} = nan(numROIs(findex),1);
    for rindex = 1:numROIs
        pre = mean(ROIs{findex}.rois(rindex).curve(pos1+control:pos2+control));
        post = mean(ROIs{findex+1}.rois(rindex).curve(pos1+control:pos2+control));
        TMI{(findex+1)/2}(rindex) = (post-pre)/(post+pre);
        ROIs{findex}.rois(rindex).TMI = TMI{(findex+1)/2}(rindex);
        ROIs{findex+1}.rois(rindex).TMI = TMI{(findex+1)/2}(rindex);
    end
end
fprintf('\tComplete.\n');


%% Save ROIdata
if saveOut && exist('ROIFiles', 'var')
    fprintf('\tSaving...\n');
    for findex = 1:numFiles
        ROIdata = ROIs{findex};
        save(ROIFiles{findex}, 'ROIdata', '-append');
        fprintf('ROIdata saved to: %s\n', ROIFiles{findex});
    end
else
    fprintf('\n');
end
