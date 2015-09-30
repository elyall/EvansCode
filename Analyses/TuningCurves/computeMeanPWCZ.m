function [PWCZmean, ROIs] = computeMeanPWCZ(ROIs, varargin)
% assumes files are input in order of [pre, post, pre, post, etc]

pos1 = 2;
pos2 = 7;
control = 1; % 1 or 0


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Start'
                pos1 = varargin{index+1};
                index = index + 2;
            case 'End'
                pos2 = varargin{index+1};
                index = index + 2;
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

%% Load data
numFiles = numel(ROIs);
fprintf('computeMeanPWCZ - %d files:', numFiles);
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


%% Compute mean of PWCZ
fprintf('\tComputing mean of PWCZ (positions %d to %d)...', pos1, pos2);
PWCZmean = cell(numFiles/2, 1);
for findex = 1:2:numFiles
    PWCZmean{(findex+1)/2} = nan(numROIs(findex),2);
    for rindex = 1:numROIs
        PWCZmean{(findex+1)/2}(rindex,1) = mean(ROIs{findex}.rois(rindex).curve(pos1+control:pos2+control));
        PWCZmean{(findex+1)/2}(rindex,2) = mean(ROIs{findex+1}.rois(rindex).curve(pos1+control:pos2+control));
        ROIs{findex}.rois(rindex).PWCZmean = PWCZmean{(findex+1)/2}(rindex,1);
        ROIs{findex+1}.rois(rindex).PWCZmean = PWCZmean{(findex+1)/2}(rindex,2);
    end
end
fprintf('\tComplete.\n');


% %% Save ROIdata
% if saveOut && exist('ROIFiles', 'var')
%     fprintf('\tSaving...\n');
%     for findex = 1:numFiles
%         ROIdata = ROIs{findex};
%         save(ROIFiles{findex}, 'ROIdata', '-append');
%         fprintf('ROIdata saved to: %s\n', ROIFiles{findex});
%     end
% else
%     fprintf('\n');
% end