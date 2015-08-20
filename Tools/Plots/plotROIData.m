function plotROIData(rois, DataInfo, ROIindex, FileIndex, FigureIndex, varargin)

saveOut = false;
saveFile = ''; % filename to save plots to

% Raster
TrialIndex = [1 inf];
frameRate = 15.45;
% showdata = true;
% showavgs = true;

% Tuning Curves
CurveColors = [0,0,1;1,0,0];
Legend = cell(2,1);
Legend{1} = [];
Legend{2} = {'Pre Control', 'Pre', 'Post Control', 'Post'};
% showFit = false;
% FitColors = [0,1,0;1,0,1];

% Constant variables
Title = [];
FigureOrder = [];


%% Parse input arguments
if ~exist('rois', 'var') || isempty(rois)
    directory = cd;
    [rois, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(rois)
        return
    elseif iscellstr(rois)
        rois = fullfile(p, rois);
    elseif ischar(rois)
        rois = {fullfile(p, rois)};
    end
end

if ~exist('ROIindex', 'var') || isempty(ROIindex)
    ROIindex = [1 inf];
end

if ~exist('FileIndex', 'var') || isempty(FileIndex)
    FileIndex = 1;
end

if ~exist('FigureIndex', 'var')
    FigureIndex = [];
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Trials','trials','TrialIndex'}
                TrialIndex = varargin{index+1};
                index = index + 2;
            case {'FrameRate', 'frameRate'}
                frameRate = varargin{index+1};
                index = index + 2;
            case 'Title'
                Title = varargin{index+1};
                index = index + 2;
            case 'FigureOrder'
                FigureOrder = varargin{index+1};
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
if ischar(rois)
    ROIFile = rois;
    if saveOut && isempty(saveFile)
        [p,fn,~] = fileparts(ROIFile);
        saveFile = fullfile(p, strcat(fn, '.pdf'));
    end
    load(ROIFile, 'ROIdata', '-mat');
    rois = ROIdata.rois;
    DataInfo = {ROIdata.DataInfo};
    clear ROIdata
elseif iscellstr(rois)
    ROIFiles = rois;
    if saveOut && isempty(saveFile)
        [p,fn,~] = fileparts(ROIFiles{1});
        saveFile = fullfile(p, strcat(fn, '.pdf'));
    end
    rois = [];
    DataInfo = cell(numel(ROIFiles), 1);
    for findex = 1:numel(ROIFiles)
        load(ROIFiles{findex}, 'ROIdata', '-mat')
        rois = cat(1, rois, ROIdata.rois);
        DataInfo{findex} = ROIdata.DataInfo;
    end
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end
totalROIs = numel(rois);


%% Determine ROIs to analyze
if ischar(ROIindex) && strcmp(ROIindex, 'all')
    ROIindex = [1, inf];
end
if ROIindex(end) == inf
    ROIindex = [ROIindex(1:end-1), ROIindex(end-1)+1:totalROIs];
end
numROIs = numel(ROIindex);

% Initialize FileIndex
if numel(FileIndex) == 1
    FileIndex = repmat(FileIndex, numROIs, 1);
end


%% Determine number of figures
if isempty(FigureIndex)
    FigureIndex = repmat(1:numROIs/2,1,2); % assumes rois = [rois(pre), rois(post)], where the rois are indexed the same
end
numFigs = max(FigureIndex);

% Initialize titles container
if isempty(Title)
    Title = cell(numFigs, 1);
end

% Determine display order
if isempty(FigureOrder)
    FigureOrder = 1:numFigs;
elseif iscolumn(FigureOrder)
    FigureOrder = FigureOrder';
end


%% Determine which trials to display
if isnumeric(TrialIndex)
    TrialIndex = repmat({TrialIndex}, numFiles);
end


%% Generate figures

% Determine y-limits for tuning curve plots
YLim = [min([rois(ROIindex).curve]), max([rois(ROIindex).curve])];

% Cycle through figures requested generating one at a time
for findex = FigureOrder
    rindices = ROIindex(FigureIndex==findex);
    
    % Determine title
    if isempty(Title{findex})
        if isempty(rois(rindices(1)).label)
            label = 'none';
        else
            label = rois(rindices(1)).label{1};
        end
        Title{findex} = sprintf('ROI: %d, Label: %s', rindices(1), label);
    end
    
    % Create figure
    hF = figure('Position', [50, 50, 1400, 800]);
    hA = zeros(1,4);
    for aindex = 1:4
        hA(aindex) = subplot(2,2,aindex);
    end
    
    % Display rasters and plot tuning curves
    CLim = zeros(2);
    for rindex = 1:2
        
        % Plot raster
        [~, CLim(rindex,:), ~] = plotIndividualRaster(...
            rois(rindices(rindex)),...
            1,...
            DataInfo{FileIndex(rindices(rindex))},...
            'axes',         hA(rindex),...
            'FrameRate',    frameRate,...
            'datatype',     'dFoF',...
            'Trials',       TrialIndex{FileIndex(rindices(rindex))});
        
        % Plot tuning curve
        plotTuningCurve(...
            rois(rindices(rindex)),...
            1,...
            'axes',             hA(4),...
            'AxesIndex',        1,...
            'curveColor',       CurveColors(rindex,:),...
            'YLim',             YLim,...
            'Legend',           Legend{rindex},...
            'Title',            Title(findex));
        
    end
    
    % Display colorbar
    axes(hA(2));
    hcb = colorbar;
    ylabel(hcb, 'dF/F');
        
    % Fix raster colormaps
    CLim = [min(CLim(:,1)), max(CLim(:,2))];
    % CLim = YLim;
    set(hA(1), 'CLim', CLim);
    set(hA(2), 'CLim', CLim);
    CMap = b2r(CLim(1), CLim(2));
    colormap(hA(1), CMap);
    colormap(hA(2), CMap);
    
    % Plot location in FoV
    axes(hA(3));
    imshow(rois(rindices(1)).mask);
    %         spatialOverlay(...
    %             rois(rindex),...
    %             [],...
    %             1,...
    %             FileIndex(index),...
    %             1,...
    %             [],...
    %             [1,1,1],...
    %             [],...
    %             'axes',         hA(FigureIndex(index),3),...
    %             'DataType',     'Discrete',...
    %             'ImgToDisplay', 'none')
    
    % Save plot to PDF
    if saveOut
        drawnow;
        export_fig(hF, saveFile, '-append');
        close(hF);
    end
    
end

