function CI = addCI(ROIdata,varargin)

N = 10000;
saveOut = false;
saveFile = '';


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'N'
                N = varargin{index+1};
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


%% Load in data
if ischar(ROIdata)
    if isempty(saveFile)
        saveFile = ROIdata;
    end
    load(ROIdata,'ROIdata','-mat');
end

Raw = gatherROIdata(ROIdata,'Raw');
[numROIs,numStim] = size(Raw);


%% Compute CIs
CI = nan(2,numStim,numROIs);
parfor_progress(numROIs);
parfor r = 1:numROIs
    for s = 1:numStim
        CI(:,s,r) = bootci(N,@mean,Raw{r,s},'type','bca');
    end
    parfor_progress;
end
parfor_progress(0);

% Distribute to ROIdata
% [ROIdata.rois(1:end).CI95] = deal(nan(2,numStim));
for r = 1:numROIs
    ROIdata.rois(r).CI95 = CI(:,:,r);
end


%% Save out
if saveOut && ~isempty(saveFile)
    if ~exist(saveFile,'file')
        save(saveFile,'ROIdata','-v7.3');
    else
        save(saveFile,'ROIdata','-append');
    end
    fprintf('Saved ROIdata to: %s\n',saveFile);
end

