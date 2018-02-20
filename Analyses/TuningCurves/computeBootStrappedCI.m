function [CI95,ROIs] = computeBootStrappedCI(Raw,varargin)

N = 10000;
ROIs = {};
ROIindex = [];
FileIndex = [];
Type = 'cper';  % 'norm', 'per', 'cper', 'bca', or 'stud' (see: BOOTCI)
Seed = 73;      % seed for random # generator
verbose = false;

saveOut = false;
saveFiles = {};


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'N'
                N = varargin{index+1};
                index = index + 2;
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
            case 'Type'
                Type = varargin{index+1};
                index = index + 2;
            case 'Seed'
                Seed = varargin{index+1};
                index = index + 2;
            case 'verbose'
                verbose = true;
                index = index + 1;
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
if ischar(saveFiles)
    saveFiles = {saveFiles};
end


%% Compute 95% CI

% Set random number generator so that outputs are reproducible
stream = RandStream('mlfg6331_64','Seed',Seed);
RandStream.setGlobalStream(stream);
% opts = statset('UseParallel',true,'UseSubstreams',true,'Streams',stream); % slower than parfor for whatever reason
opts = statset('Streams',stream);

[numROIs,numStim] = size(Raw);
fprintf('Computing bootstrapped confidence intervals for %d ROIs...',numROIs);

CI95 = nan(2,numStim,numROIs);
if verbose
    p = parfor_progress(numROIs);
end
parfor ind = 1:numROIs
    for sind = 1:numStim
        reset(stream); % reset random # generator
        CI95(:,sind,ind) = bootci(N,{@mean,Raw{ind,sind}},'type',Type,'Options',opts); % bootstrapped confidence intervals of the mean
    end
    if verbose
        parfor_progress(p);
    end
end
if verbose
    parfor_progress(p,0);
end

fprintf('\tComplete\n');


%% Distribute to ROIdata structure
if (nargout>1 || (saveOut && ~isempty(saveFiles))) && ~isempty(ROIs)
    for ind = 1:numROIs
        ROIs{FileIndex(ind)}.rois(ROIindex(ind)).CI95 = CI95(:,:,ind);
    end
    if saveOut
        for ind = 1:numel(ROIs)
            ROIdata = ROIs{ind};
            fprintf('Saving ROIdata to: %s...',saveFiles{ind});
            if ~exist(saveFiles{ind},'file')
                save(saveFiles{ind},'ROIdata','-v7.3');
            else
                save(saveFiles{ind},'ROIdata','-append');
            end
            fprintf('\tComplete\n');
        end
    end
end

CI95 = permute(CI95,[3,1,2]); % index ROIs in first dimension