function [PValues, ROIdata] = tuningAnova(ROIdata, varargin)


StimIndex = [2 inf];
ROIindex = [1 inf];

saveOut = false;
saveFile = '';


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
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

if ~exist('ROIdata', 'var') || isempty(ROIdata)
    [ROIdata,p] = uigetfile({'*.rois;*.mat'}, 'Select ROI files:', directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = {fullfile(p, ROIdata)};
end


%% Load data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIFile, 'ROIdata', '-mat');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end


%% Determine ROIs to analyze
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:numel(ROIdata.rois)];
end
numROIs = numel(ROIindex);


%% Determine stimuli to analyze
numStims = length(ROIdata.rois(ROIindex(1)).curve);
if StimIndex(end) == inf
    StimIndex = [StimIndex(1:end-1), StimIndex(end-1)+1:numStims];
end
numStims = numel(StimIndex);


%% Compute anova

% Initialize output
PValues = nan(numROIs, 1); % [btwn positions, btwn pre and post, interaction btwn the two]

% Compute significance
parfor_progress(numROIs);
for rindex = 1:numROIs
    
    % Gather data
    temp = ROIdata.rois(ROIindex(rindex)).Raw;
    N = cellfun(@numel, temp);
    N = max(N(:));
    
    % Organize data
    dFoF = [];
    for tindex = StimIndex
        dFoF = cat(1, dFoF, [temp{tindex}; nan(N-numel(temp{tindex}),1)]);
    end
    
    % Create dictionary
    g1 = reshape(repmat(1:numStims, N, 1), numStims*N, 1);
    
    % Compute 2-way ANOVA
    PValues(rindex) = anovan(dFoF, g1, 'model','full','display','off');
    
    parfor_progress;
end
parfor_progress(0);


%% Distribute to struct
if nargin > 1 || saveOut
    for rindex = 1:numROIs
        ROIdata.rois(ROIindex(rindex)).TunedPValue = PValues(rindex);
    end
end


%% Save to file
if saveOut
    if ~exist(saveFile, 'file')
        save(saveFile, 'ROIdata', '-mat', '-v7.3');
    else
        save(saveFile, 'ROIdata', '-mat', '-append');
    end
    fprintf('\tROIdata saved to: %s\n', saveFile);
end

