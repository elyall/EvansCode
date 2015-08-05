function rois = plotROIData(rois, DataInfo, FileIndex, ROIindex, varargin)
% ROIFiles are input {pre1, post1, pre2, post2, ...}

saveOut = false;
saveFile = ''; % filename to save plots to

% Raster
TrialIndex = [1 inf];
showdata = true;
showavgs = true;
frameRate = 15.45;

% Tuning Curves
showFit = false;
curveColor1 = 'r';
fitColor1 = 'g';
curveColor2 = 'b';
fitColor2 = 'm';


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

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Trials','trials','TrialIndex'}
                TrialIndex = varargin{index+1};
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
elseif iscellstr(rois)
    ROIFiles = rois;
    if saveOut && isempty(saveFile)
        [p,fn,~] = fileparts(ROIFiles{1});
        saveFile = fullfile(p, strcat(fn, '.pdf'));
    end
    rois = [];
    for findex = 1:numel(ROIFiles)
        load(ROIFiles{findex}, 'ROIdata', '-mat')
        rois = cat(1, rois, ROIdata.rois);
    end
end
if saveOut && isempty(saveFile)
    warning('Cannot save output as no file specified');
    saveOut = false;
end


%% Determine matched ROIs (doesn't account for double counted ROIs in border region)
gi


%% Order ROIs by distance from point selected
Order = orderROIsByLocation(rois(1:2:end), DataInfo(1:2:end), ROIindex(:,1), (FileIndex(:,1)+1)/2);
ROIindex = ROIindex(Order, :);
FileIndex = FileIndex(Order, :);


%% Determine which trials to display
if isnumeric(TrialIndex)
    TrialIndex = repmat({TrialIndex}, numFiles);
end


%% Generate figures
YLim = [min([rois(:).curve]), max([rois(:).curve])];

% Cycle through ROIs in file generate figures one at a time
for rindex = 1:numROIs
    
    % Create figure
    hA = cell(1, 4);
    hF = figure('Position', [50, 50, 1400, 800]);
    for aindex = 1:4
        hA{1,aindex} = subplot(2,2,aindex);
    end
    
    % raster
    [~, CLim, ~] = plotIndividualRaster(...
        rois{FileIndex(rindex,1)}.rois(ROIindex(rindex,1)),...
        1,...
        rois{FileIndex(rindex,1)}.DataInfo,...
        'axes', hA(:,1),...
        'frameRate', frameRate,...
        'Trials', TrialIndex{FileIndex(rindex,1)});
    
    % tuning curve
    plotTuningCurve2(...
        rois{FileIndex(rindex,1)}.rois(ROIindex(rindex,1)),...
        1,...
        'axes', hA(:,4),...
        'fitColor', fitColor1,...
        'curveColor', curveColor1,...
        'showFit',...
        'legendLocation', 'NorthWest');
    
    % LOCATION
    axes(hA{1, 3});
    spatialOverlay(...
        rois,...
        DataInfo,...
        ROIindex(rindex,1),...
        FileIndex(rindex,1),...
        1,...
        [],...
        [1,1,1],...
        [],...
        'axes', hA{1,3},...
        'DataType', 'Discrete',...
        'ImgToDisplay', 'none')
    %         imshow(ROIs{FileIndex(rindex,1)}.rois(ROIindex(rindex,1)).mask);
    
    % POST-TRIMMING
    % raster
    [~, CLim, cmap] = plotIndividualRaster(...
        rois{FileIndex(rindex,2)}.rois(ROIindex(rindex,2)),...
        1,...
        rois{FileIndex(rindex,1)}.DataInfo,...
        'axes', hA(:,2),...
        'frameRate', 15.45,...
        'Trials', TrialIndex{FileIndex(rindex,1)},...
        'clim', CLim,...
        'showColorBar');
    
    % tuning curve
    plotTuningCurve(...
        rois{FileIndex(rindex,2)}.rois(ROIindex(rindex,2)),...
        1,...
        'axes', hA(:,4),...
        'fitColor', fitColor2,...
        'curveColor', curveColor2,...
        'showFit',...
        'legendLocation', 'NorthEast',...
        'Title', {sprintf('ROI %s (%d)',...
        Tags{rindex}, ROIindex(rindex))});
    
    % Fix pre-trimming colormaps
    set(hA{1, 1}, 'CLim', CLim);
    colormap(hA{1, 1}, cmap{1});
    
    if saveOut
        drawnow
        export_fig(saveFile, hF, '-append');
        close(hF);
    end
    
end

