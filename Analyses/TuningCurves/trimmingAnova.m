function PValues = trimmingAnova(ROIs, ROIindex, varargin)
% ROIs - 1x2 or 2x1 cell array of ROIdata structs or cell array of filenames
% ROIindex - Nx2 array of ROI indices to analyze or cell array of filenames
% to load Maps from for auto matching or vector of imref2d w/ numel=2
%
% PValues - Nx3 array of p-values -> [btwn positions, btwn pre and post, interaction btwn the two]

StimIndex = [2 inf];


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'StimIndex'
                StimIndex = varargin{index+1};
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

if ~exist('ROIs', 'var') || isempty(ROIs)
    [ROIs,p] = uigetfile({'*.rois;*.mat'}, 'Select ROI files:', directory, 'MultiSelect', 'on');
    if isnumeric(ROIs)
        return
    elseif iscellstr(ROIs)
        ROIs = fullfile(p, ROIs);
    else
        ROIs = {fullfile(p, ROIs)};
    end
end


%% Load data
if iscellstr(ROIs)
    ROIFiles = ROIs;
    for findex = 1:numel(ROIFiles)
        load(ROIFiles{1}, 'ROIdata', '-mat');
        ROIs{findex} = ROIdata;
    end
end

% Determine ROIs to analyze
if ~exist('ROIindex', 'var') || isempty(ROIindex)
    numROIs = cellfun(@(x) (numel(x.rois)), ROIs); % vector of number of rois in each element of cell array
    if numel(unique(numROIs))>1
        error('Both files must have the same number of ROIs');
    end
    ROIindex = repmat((1:numROIs(1))',1,2);
elseif iscellstr(ROIindex) || isa(ROIindex, 'imref2d')
    ROIindex = autoMatchROIs(ROIs, ROIindex);
end
numROIs = size(ROIindex, 1);


%% Determine stimuli to analyze
if StimIndex(end) == inf
    StimIndex = [StimIndex(1:end-1), StimIndex(end-1)+1:numel(ROIs{1}.rois(ROIindex(1,1)).curve)];
end
numStims = numel(StimIndex);


%% Compute anova

% Initialize output
PValues = nan(numROIs, 3); % [btwn positions, btwn pre and post, interaction btwn the two]

% Compute significance
parfor_progress(numROIs);
for rindex = 1:numROIs
    
    % Gather data
    temp1 = ROIs{1}.rois(ROIindex(rindex,1)).Raw;
    temp2 = ROIs{2}.rois(ROIindex(rindex,2)).Raw;
    N = cellfun(@numel, [temp1,temp2]);
    N = max(N(:));
    
    % Organize data
    dFoFPre = [];
    dFoFPost = [];
    for tindex = StimIndex
        dFoFPre = cat(1, dFoFPre, [temp1{tindex}; nan(N-numel(temp1{tindex}),1)]);
        dFoFPost = cat(1, dFoFPost, [temp2{tindex}; nan(N-numel(temp2{tindex}),1)]);
    end
    
    % Create dictionary
    g1 = repmat(reshape(repmat(1:numStims, N, 1), numStims*N, 1), 2, 1);
    g2 = reshape(repmat([numStims+1, numStims+2], numStims*N, 1), numStims*N*2, 1);
    
    % Compute 2-way ANOVA
    PValues(rindex,:) = anovan([dFoFPre;dFoFPost], [g1,g2], 'model','full','display','off')';
    
    parfor_progress;
end
parfor_progress(0);

