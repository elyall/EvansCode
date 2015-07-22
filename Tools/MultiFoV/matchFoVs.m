function [Maps, ExperimentFiles] = matchFoVs(ExperimentFiles)

%% Check arguments
dataAvailable = false;
if exist('ExperimentFiles', 'var') && ~isempty(ExperimentFiles)
    dataAvailable = true;
end

%% Configure default settings
gd.Internal.directory = 'D:\Evan\Data';
gd.Internal.Settings.ROITextSize = 20;
gd.Internal.Settings.ROILineWidth = 2;

%% Create & Populate Figure

% FIGURE
if nargout
    gd.fig = figure(...
        'NumberTitle',          'off',...
        'Name',                 'Match Fields of View',...
        'ToolBar',              'none',...
        'Units',                'pixels',...
        'Position',             [50, 50, 1400, 800],...
        'KeyPressFcn',          @(hObject,eventdata)KeyPressCallback(hObject,eventdata,guidata(hObject)),...
        'CloseRequestFcn',      'uiresume(gcbf)');
else
    gd.fig = figure(...
        'NumberTitle',          'off',...
        'Name',                 'Match Fields of View',...
        'ToolBar',              'none',...
        'Units',                'pixels',...
        'Position',             [50, 50, 1400, 800],...
        'KeyPressFcn',          @(hObject,eventdata)KeyPressCallback(hObject,eventdata,guidata(hObject)));
end

% DISPLAY
% panel
gd.Display.panel = uipanel(...
    'Title',                'Merged Display',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [0,0,.7,1]);
% axes
gd.Display.axes = axes(...
    'Parent',               gd.Display.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,1,1]);
axis off

% CONTROLS
% panel
gd.Control.panel = uipanel(...
    'Title',                'Controls',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [.7,0,.3,1]);
% load button
gd.Control.load = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Load Data',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,.9,.5,.1],...
    'Callback',             @(hObject,eventdata)LoadExperiment(hObject,eventdata,guidata(hObject)),...
    'BackgroundColor',      [.3,.3,.3],...
    'ForegroundColor',      [1,1,1]);
% save button
gd.Control.save = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Save Offsets',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [.5,.9,.5,.1],...
    'Callback',             @(hObject,eventdata)SaveData(hObject,eventdata,guidata(hObject)),...
    'BackgroundColor',      [.7,1,.7],...
    'Enable',               'off');
% popupmenu: merge type
gd.Control.overlay = uicontrol(...
    'Style',                'popupmenu',...
    'String',               {'Red/Green','Green/Magenta','Blend','Difference'},...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,.87,1,.03],...
    'Callback',             @(hObject,eventdata)plotDataAxes(guidata(hObject)));
% checkbox: invert
gd.Control.invert = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Invert?',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,.85,.25,.02],...
    'Callback',             @(hObject,eventdata)plotDataAxes(guidata(hObject)));
% table: cnames= name, view (logical), move (logical), side (popup: A or B)
gd.Control.selection = uitable(...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,1,.85],...
    'RowName',              [],...
    'ColumnName',           {'File','View','Move','Y','X','Path','Remove'},...
    'ColumnEditable',       [false,true,true,false,false,false,true],...
    'ColumnFormat',         {'char','logical','logical','char','char','char','logical'},...
    'ColumnWidth',          {100,40,40,40,40,200,40},...
    'Enable',               'off',...
    'CellEditCallback',     @(hObject,eventdata)DataSelection(hObject,eventdata,guidata(hObject)));

guidata(gd.fig, gd);

% Load files
if dataAvailable
    LoadExperiment(gd.Control.load,ExperimentFiles,gd)
end

% Return requested outputs
if nargout
    uiwait(gd.fig);
    gd = guidata(gd.fig);
    if isfield(gd, 'Experiment')
        Maps = {gd.Experiment(:).Map};
        ExperimentFiles = {gd.Experiment(:).filename};
    else
        Maps = {};
        ExperimentFiles = {};
    end
    delete(gd.fig)
end

%% KEYPRESS CALLBACKS
function KeyPressCallback(hObject,eventData,gd)

% Determine step size
if isempty(eventData.Modifier)
    step = 1;
elseif numel(eventData.Modifier)==1
    step = 10;
else
    step = 100;
end

% Select images
contents = get(gd.Control.selection, 'Data');
viewSelection = find([contents{:,2}]);
indices = viewSelection([contents{viewSelection, 3}]);

% Move images
switch eventData.Key
    case {'uparrow','w'} %translate ROI up
        for index = indices
            gd.Experiment(index).Map.YWorldLimits = gd.Experiment(index).Map.YWorldLimits - step;
            contents{index, 4} = gd.Experiment(index).Map.YWorldLimits(1) - 0.5;
        end
    case {'downarrow','s'} %translate ROI down
        for index = indices
            gd.Experiment(index).Map.YWorldLimits = gd.Experiment(index).Map.YWorldLimits + step;
            contents{index, 4} = gd.Experiment(index).Map.YWorldLimits(1) - 0.5;
        end
    case {'rightarrow','d'} %translate ROI right
        for index = indices
            gd.Experiment(index).Map.XWorldLimits = gd.Experiment(index).Map.XWorldLimits + step;
            contents{index, 5} = gd.Experiment(index).Map.XWorldLimits(1) - 0.5;
        end
    case {'leftarrow','a'} %translate ROI left
        for index = indices
            gd.Experiment(index).Map.XWorldLimits = gd.Experiment(index).Map.XWorldLimits - step;
            contents{index, 5} = gd.Experiment(index).Map.XWorldLimits(1) - 0.5;
        end
end
guidata(hObject, gd);
set(gd.Control.selection, 'Data', contents);
plotDataAxes(gd);

%% MAIN
function LoadExperiment(hObject,eventdata,gd)

% Select Files
if ~iscell(eventdata)
    if isdir(gd.Internal.directory)
        directory = gd.Internal.directory;
    elseif exist('CanalSettings.m', 'file');
        directory = CanalSettings('ExperimentDirectory');
    else
        directory = cd;
    end
    ExperimentFiles = uipickfiles('Prompt', 'Choose experiment files', 'FilterSpec', [directory, '*.ciexp']);
    if isnumeric(ExperimentFiles)
        return
    end
else
    ExperimentFiles = eventdata;
end
[gd.Internal.directory,~,~] = fileparts(ExperimentFiles{1});
numFiles = numel(ExperimentFiles);

% Load Experiment Data
if isfield(gd, 'Experiment') % previous files loaded
    offset = numel(gd.Experiment); % add files to end
else
    offset = 0;
end
for findex = 1:numFiles
    gd.Experiment(offset + findex).filename = ExperimentFiles{findex};
    load(ExperimentFiles{findex}, 'Map', 'ImageFiles', '-mat');
    if exist('ImageFiles', 'var')
        gd.Experiment(offset + findex).ImageFiles = ImageFiles;
    else
        error('Compute average of data first using ''ComputeProjections''');
    end
    if exist('Map', 'var')
        gd.Experiment(offset + findex).Map = Map;
    else
        gd.Experiment(offset + findex).Map = [];
    end
    
end

% Build Individual Images
for findex = offset+1:offset+numFiles
    gd.Experiment(findex).Image = permute(gd.Experiment(findex).ImageFiles.Average, [1,2,4,3]);
    [H, W, C, ~] = size(gd.Experiment(findex).Image);
    if C == 1
        gd.Experiment(findex).Image = cat(3, zeros(H, W), gd.Experiment(findex).Image, zeros(H, W)); % green only
    elseif C == 2
        gd.Experiment(findex).Image = cat(3, gd.Experiment(findex).Image(:,:,[2,1]), zeros(H, W)); % green & red data
    end
    if isempty(gd.Experiment(findex).Map)
        gd.Experiment(findex).Map = imref2d([H,W]);
    end
    % Determine colormaps
    for cindex = 2
        temp = gd.Experiment(findex).Image(:,:,cindex);
        gd.Experiment(findex).clim(cindex).limits = [min(temp(:)), max(temp(:))];
        gd.Experiment(findex).clim(cindex).current = [gd.Experiment(findex).clim(cindex).limits(1), prctile(temp(:),99.98)];
        set(gd.Display.axes,'CLim',gd.Experiment(findex).clim(cindex).current);
        gd.Experiment(findex).colormap(cindex).current = colormap;
    end
end

% Set current colormaps
gd = AdjustColorMaps(gd, offset+1:offset+numFiles);

% Add Files to Table
contents = get(gd.Control.selection, 'Data');
for findex = offset+1:offset+numFiles %{'File','View','Move','Group','X','Y','Path'}
    [p,f,~] = fileparts(gd.Experiment(findex).filename);
    contents{findex, 1} = f;
    contents{findex, 2} = false;
    contents{findex, 3} = false;
    contents{findex, 4} = gd.Experiment(findex).Map.YWorldLimits(1) - 0.5;
    contents{findex, 5} = gd.Experiment(findex).Map.XWorldLimits(1) - 0.5;
    contents{findex, 6} = p;
end
contents{1, 2} = true; % view first file
set(gd.Control.selection, 'Data', contents, 'Enable', 'on');
set(gd.Control.save, 'Enable', 'on');

guidata(hObject, gd);
plotDataAxes(gd);

function gd = AdjustColorMaps(gd, indices)
for findex = indices
    gd.Experiment(findex).current = gd.Experiment(findex).Image;
    % Adjust colormap
    for cindex = 2
        temp = gd.Experiment(findex).current(:,:,cindex);
        temp = temp./gd.Experiment(findex).clim(cindex).current(2);
        temp(temp>1) = 1;
        temp(temp<gd.Experiment(findex).clim(cindex).current(1)/gd.Experiment(findex).clim(cindex).current(2)) = gd.Experiment(findex).clim(cindex).current(1)/gd.Experiment(findex).clim(cindex).current(2);
        gd.Experiment(findex).current(:,:,cindex) = temp;
    end
end

function gd = plotDataAxes(gd)

% Select images
contents = get(gd.Control.selection, 'Data');
viewSelection = find([contents{:,2}]);
fixed = find(~[contents{viewSelection, 3}]);
moving = find([contents{viewSelection, 3}]);

% Determine if colormap should be inverted
if get(gd.Control.invert,'Value')
    invert = true;
else
    invert = false;
end

% Create two distinct images
if ~isempty(fixed)
    if ~invert
        ImageA = gd.Experiment(viewSelection(fixed(1))).current(:,:,2);
    else
        ImageA = 1 - gd.Experiment(viewSelection(fixed(1))).current(:,:,2);
    end
    MapA = gd.Experiment(viewSelection(fixed(1))).Map;
    for index = 2:numel(fixed)
        if ~invert
            ImageA2 = gd.Experiment(viewSelection(fixed(index))).current(:,:,2);
        else
            ImageA2 = 1 - gd.Experiment(viewSelection(fixed(index))).current(:,:,2);
        end
        [ImageA, MapA] = imfuse(...
            ImageA,...
            MapA,...
            ImageA2,...
            gd.Experiment(viewSelection(fixed(index))).Map,...
            'blend',...
            'Scaling', 'Independent');
    end
else
    ImageA = [];
end
if ~isempty(moving)
    if ~invert
        ImageB = gd.Experiment(viewSelection(moving(1))).current(:,:,2);
    else
        ImageB = 1 - gd.Experiment(viewSelection(moving(1))).current(:,:,2);
    end
    MapB = gd.Experiment(viewSelection(moving(1))).Map;
    for index = 2:numel(moving)
        if ~invert
            ImageB2 = gd.Experiment(viewSelection(moving(index))).current(:,:,2);
        else
            ImageB2 = 1 - gd.Experiment(viewSelection(moving(index))).current(:,:,2);
        end
        [ImageB, MapB] = imfuse(...
            ImageB,...
            MapB,...
            ImageB2,...
            gd.Experiment(viewSelection(moving(index))).Map,...
            'blend',...
            'Scaling', 'Independent');
    end
else
    ImageB = [];
end

% Determine display type
col = [1,1,1];
DisplayType = get(gd.Control.overlay, 'Value');
if DisplayType == 1
    DisplayType = 'falsecolor';
    col = [2,1,2];
elseif DisplayType == 2
    DisplayType = 'falsecolor';
    col = [1,2,2];
elseif DisplayType == 3
    DisplayType = 'blend';
elseif DisplayType == 4
    DisplayType = 'diff';
end

% Create merged image
if ~isempty(ImageA) && ~isempty(ImageB)
    switch DisplayType
        case 'falsecolor'
            img = imfuse(ImageA,MapA,ImageB,MapB,DisplayType,'Scaling','Independent','ColorChannels',col);
        otherwise
            img = imfuse(ImageA,MapA,ImageB,MapB,DisplayType,'Scaling','Independent');
    end
elseif isempty(ImageB)
    img = ImageA;
else
    img = ImageB;
end

% Select axes
axes(gd.Display.axes);

% Display Image
imghandle = imagesc(img);
% colormap([gd.Experiment(dataIndex).colormap(1).current(:,1),...
%     gd.Experiment(dataIndex).colormap(2).current(:,2),...
%     (0:1/(size(gd.Experiment(dataIndex).colormap(1).current,1)-1):1)']);
set(gca,'xtick',[],'ytick',[])

guidata(gd.fig,gd); % update guidata

function DataSelection(hObject,eventdata,gd)

if eventdata.Indices(2) == 2 || eventdata.Indices(2) == 3 % display selected images
    plotDataAxes(gd);
    
elseif eventdata.Indices(2) == 7 % remove selected file
    gd.Experiment(eventdata.Indices(1)) = [];
    contents = get(hObject, 'Data');
    contents(eventdata.Indices(1),:) = [];
    set(hObject, 'Data', contents);
    guidata(hObject, gd);
    plotDataAxes(gd);
end
% error handling for input data

function SaveData(hObject,eventdata,gd)
answer = questdlg('Are you sure? Any past Maps will be overwritten');
if ~strcmp(answer, 'Yes')
    return
end
set(gcf, 'pointer', 'watch')
drawnow
for findex = 1:numel(gd.Experiment);
    Map = gd.Experiment(findex).Map;
    save(gd.Experiment(findex).filename, 'Map', '-append');
    fprintf('Map saved to: %s\n', gd.Experiment(findex).filename)
end
set(gcf, 'pointer', 'arrow')

% function CloseReqFcn(src, ~)
% pause(.5);
% delete(gcf);
