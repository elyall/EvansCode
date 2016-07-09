function plotROIChange(ROIs, ROIindex, FileIndex, varargin)

% Saving
saveOut = false;
saveFile = '';                  % filename to save plots to

% Tuning Curves
Title = {};
PWCZ = [];

% Raster
TrialIndex = [];
stimorder = [];

% Figure properties
y = 4;                          % # of rois per fig
FigPos = [100,100,1200,800];    % position in pixels


%% Parse input arguments
if ~exist('ROIs', 'var') || isempty(ROIs)
    directory = cd;
    [ROIs, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(ROIs)
        return
    else
        ROIs = fullfile(p, ROIs);
    end
end

if ~exist('ROIindex', 'var') || isempty(ROIindex)
    ROIindex = [];
end
if ~exist('FileIndex', 'var') || isempty(FileIndex)
    FileIndex = [];
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Title'
                Title = varargin{index+1};
                index = index + 2;
            case 'PWCZ'
                PWCZ = varargin{index+1};
                index = index + 2;
            case {'Trials','trials','TrialIndex'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case 'stimorder'
                stimorder = varargin{index+1};
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
if iscellstr(ROIs)
    ROIFiles = ROIs;
    for f = 1:numel(ROIFiles)
        load(ROIFiles{f}, 'ROIdata', '-mat')
        ROIs{f} = ROIdata;
    end
    if saveOut && isempty(saveFile)
        [p,fn,~] = fileparts(ROIFiles{1});
        saveFile = fullfile(p,[fn,'.pdf']);
    end
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end
numFiles = numel(ROIs);


%% Determine ROIs to analyze
if isempty(ROIindex)
    FileIndex = [];
    for f=1:2:numFiles
        Min = min(cellfun(@(x) numel(x.rois),ROIs(f:f+1)));
        ROIindex = [ROIindex;repmat((1:Min)',1,2)];
        FileIndex = [FileIndex;repmat([f,f+1],Min,1)];
    end
elseif ~isequal(size(ROIindex),size(FileIndex))
    error('ROIindex and FileIndex have to be the same size');
end
numROIs = size(ROIindex,1);

% Set titles
if isempty(Title)
    Title = cell(numROIs,1);
end


%% Determine which trials to display
if isempty(TrialIndex)
    TrialIndex = cell(numFiles,1);
    for f = 1:numFiles
        TrialIndex{f} = ROIs{f}.DataInfo.TrialIndex;
    end
end


%% Display data
x=3;
for r = 1:numROIs
    
    % Create figure
    if mod(r,y)==1
        figure('Position',FigPos);
    end
    
    % Subplot index
    s = mod(x*(r-1)+1,x*y);
    
    %% Tuning Curves
    
    % Determine Limits
    YLim = nan(1,2);
    YLim(1) = min([ROIs{FilesIndex(r,1)}.rois(ROIindex(r,1)).curve(2:end),ROIs{FilesIndex(r,2)}.rois(ROIindex(r,2)).curve(2:end)]);
    YLim(2) = max([ROIs{FilesIndex(r,1)}.rois(ROIindex(r,1)).curve(2:end),ROIs{FilesIndex(r,2)}.rois(ROIindex(r,2)).curve(2:end)]);
    YLim = 1.2*YLim;
    
    % Plot curves
    plotTuningCurve(ROIs{FilesIndex(r,1)}.rois, ROIindex(r,1),...
        'curveColor', 'k', 'Title', '', 'axes', subplot(y,x,s), 'AxesIndex', 1);
    plotTuningCurve(ROIs{FilesIndex(r,2)}.rois, ROIindex(r,2),...
        'curveColor', 'r', 'Title', '', 'axes', subplot(y,x,s), 'AxesIndex', 1,...
        'PWCZ', PWCZ, 'YLim', YLim, 'Title', Title{r});
    
    %% Rasters
    CLim = zeros(2);
    [~, CLim(1,:), ~] = plotIndividualRaster(ROIs{FilesIndex(r,1)}.rois,ROIindex(r,1),...
        ROIs{FilesIndex(r,1)}.DataInfo,...
        'axes',         subplot(y,x,s+1),...
        'datatype',     'dFoF',...
        'Trials',       TrialIndex{1},...
        'stimOrder',    stimorder);
    [~, CLim(2,:), ~] = plotIndividualRaster(ROIs{FilesIndex(r,2)}.rois,ROIindex(r,2),...
        ROIs{FilesIndex(r,2)}.DataInfo,...
        'axes',         subplot(y,x,s+2),...
        'datatype',     'dFoF',...
        'Trials',       TrialIndex{2},...
        'stimOrder',    stimorder);
    CLim = [min(CLim(:,1)), max(CLim(:,2))];
    set(subplot(y,x,s+1), 'CLim', CLim);
    set(subplot(y,x,s+2), 'CLim', CLim);
    CMap = b2r(CLim(1), CLim(2));
    colormap(subplot(y,x,s+1), CMap);
    colormap(subplot(y,x,s+2), CMap);
    
    %% Export
    if saveOut && ~isempty(saveFile) && (~mod(r,y) || r == numROIs)
        export_fig(gcf, saveFile, '-append');
        close(gcf);
    end
    
end

