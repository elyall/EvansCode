function [CI95,ROIs] = computeBootStrappedCI(Raw,varargin)

ROIs = {};
ROIindex = [];
FileIndex = [];

saveOut = false;
saveFiles = {};


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIs'
                ROIs = varargin{index+1};
                index = index + 2;
                if ~iscell(ROIs)
                    ROIs = {ROI};
                end
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FileIndex'
                FileIndex = varargin{index+1};
                index = index + 2;
            case {'save','Save'}
                saveOut = true;
                index = index + 1;
            case 'saveFiles'
                saveFiles = varargin{index+1};
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

if ischar(Raw) || iscellstr(Raw)
    ROIFiles = Raw;
    if saveOut && isempty(saveFiles)
        saveFiles = ROIFiles;
    end
    [Raw,FileIndex,ROIindex,~,~,ROIs] = gatherROIdata(ROIFiles,'Raw');
end


%% Compute 95% CI

[numROIs,numStim] = size(Raw);
fprintf('Computing bootstrapped confidence intervals for %d ROIs...',numROIs);

CI95 = nan(numROIs,2,numStim);
parfor ind = 1:numROIs
    for sind = 1:numStim
        CI95(ind,:,sind) = bootci(10000,{@mean,Raw{ind,sind}},'type','bca'); % bootstrapped confidence intervals of the mean
    end
%     b = cellfun(@(x) bootci(10000,{@mean,x},'type','bca'), Raw(ind,:),'UniformOutput',false);
%     CI95(ind,:,:) = cat(2,b{:});
end

fprintf('\tComplete\n');


%% Distribute to ROIdata structure
if (nargout>1 || (saveOut && ~isempty(saveFiles))) && ~isempty(ROIs)
    for ind = 1:numROIs
        ROIs{FileIndex(ind)}.rois(ROIindex(ind)).CI95 = squeeze(CI95(ind,:,:));
    end
    if saveOut
        for ind = 1:numel(ROIs)
            ROIdata = ROIs{ind};
            if ~exist(saveFiles{ind},'file')
                save(saveFiles{ind},'ROIdata','-v7.3');
            else
                save(saveFiles{ind},'ROIdata','-append');
            end
        end
    end
end