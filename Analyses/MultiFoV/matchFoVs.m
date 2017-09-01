function [Maps, Images] = matchFoVs(Images, Maps)
% Images is a cell array of N images or a matrix of HxWxN
% Maps is a structure of imref2d registration files


gd.Internal.directory = cd;

%% Parse input arguments
dataAvailable = false;
if exist('Images', 'var') && ~isempty(Images)
    dataAvailable = true;
    
    if isnumeric(Images)
        temp = Images;
        Images = cell(size(temp,3),1);
        for iindex = 1:size(temp,3)
            Images{iindex} = temp(:,:,iindex);
        end
    end
end


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
    'ColumnEditable',       [true,true,true,true,true,true,true],...
    'ColumnFormat',         {'char','logical','logical','char','char','char','logical'},...
    'ColumnWidth',          {100,40,40,40,40,200,40},...
    'Enable',               'off',...
    'CellEditCallback',     @(hObject,eventdata)DataSelection(hObject,eventdata,guidata(hObject)));

guidata(gd.fig, gd);

% Load files
if dataAvailable
    if ~exist('Maps', 'var') || isempty(Maps)
        Maps = [];
    end
    LoadExperiment(gd.Control.load,{Images, Maps},gd)
end

% Return requested outputs
if nargout
    uiwait(gd.fig);
    gd = guidata(gd.fig);
    if isfield(gd, 'DataSets')
        temp = {gd.DataSets(:).Map};
        Maps = imref2d();
        for f = 1:numel(temp)
            Maps(f) = temp{f};
        end
        Images = {gd.DataSets(:).filename};
    else
        Maps = {};
        Images = {};
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
            gd.DataSets(index).Map.YWorldLimits = gd.DataSets(index).Map.YWorldLimits - step;
            contents{index, 4} = gd.DataSets(index).Map.YWorldLimits(1) - 0.5;
        end
    case {'downarrow','s'} %translate ROI down
        for index = indices
            gd.DataSets(index).Map.YWorldLimits = gd.DataSets(index).Map.YWorldLimits + step;
            contents{index, 4} = gd.DataSets(index).Map.YWorldLimits(1) - 0.5;
        end
    case {'rightarrow','d'} %translate ROI right
        for index = indices
            gd.DataSets(index).Map.XWorldLimits = gd.DataSets(index).Map.XWorldLimits + step;
            contents{index, 5} = gd.DataSets(index).Map.XWorldLimits(1) - 0.5;
        end
    case {'leftarrow','a'} %translate ROI left
        for index = indices
            gd.DataSets(index).Map.XWorldLimits = gd.DataSets(index).Map.XWorldLimits - step;
            contents{index, 5} = gd.DataSets(index).Map.XWorldLimits(1) - 0.5;
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
    else
        directory = cd;
    end
    [Images, p] = uigetfile({'*.exp;*.align'}, 'Select files to register:', directory, 'MultiSelect', 'on');
    if isnumeric(Images)
        return
    end
    gd.Internal.directory = p;
    guidata(hObject, gd);
    if iscellstr(Images)
        Images = fullfile(p, Images);
    elseif ischar(Images)
        Images = {fullfile(p, Images)};
    end
    Maps = cell(numel(Images), 1);
else
    Images = eventdata{1};
    if ~isempty(eventdata{2})
        Maps = eventdata{2};
    else
        Maps = cell(numel(Images), 1);
        for findex = 1:numel(Images)
            Maps{findex} = imref2d([size(Images{findex},1), size(Images{findex},2)]);
        end
    end
end

% Load Images
if iscellstr(Images)
    ImageFiles = Images;
    Images = cell(numel(ImageFiles),1);
    for findex = 1:numel(ImageFiles)
        [~,~,ext] = fileparts(ImageFiles{findex});
        switch ext
            case '.align'
                temp = load(ImageFiles{findex}, 'm', 'Map', '-mat');
                Images{findex} = temp.m;
            case '.exp'
                temp = load(ImageFiles{findex}, 'ImageFiles', 'Map', '-mat');
                Images{findex} = squeeze(temp.ImageFiles.Average(:,:,1,1));
        end
        if isfield(temp, 'Map')
            Maps{findex} = temp.Map;
        else
            Maps{findex} = imref2d([size(Images{findex},1), size(Images{findex},2)]);
        end
    end
else
    ImageFiles = cellfun(@num2str, num2cell(1:numel(Images)), 'UniformOutput', false);
    ImageFiles = fullfile(gd.Internal.directory, strcat('file', ImageFiles, '.exp'));
end
numFiles = numel(Images);

% Initialize struct
if ~isfield(gd, 'DataSets') % previous files loaded
    gd.DataSets = [];
    offset = 0;
else
    offset = numel(gd.DataSets);
end

% Add data to struct
for findex = 1:numFiles
    
    gd.DataSets(offset+findex).filename = ImageFiles{findex};
    gd.DataSets(offset+findex).Image = Images{findex};
    gd.DataSets(offset+findex).Map = Maps{findex};

    % Determine colormap
    gd.DataSets(offset+findex).clim.limits = [min(gd.DataSets(offset+findex).Image(:)), max(gd.DataSets(offset+findex).Image(:))];
    gd.DataSets(offset+findex).clim.current = [gd.DataSets(offset+findex).clim.limits(1), prctile(gd.DataSets(offset+findex).Image(:),99.98)];
    set(gd.Display.axes, 'CLim', gd.DataSets(offset+findex).clim.current);
    gd.DataSets(offset+findex).colormap.current = colormap;
    
    % Create image
    gd.DataSets(offset+findex).display = gd.DataSets(offset+findex).Image;
end

% Set current colormaps
% gd = AdjustColorMaps(gd, offset+1:offset+numFiles);

% Add Files to Table
contents = get(gd.Control.selection, 'Data');
for findex = offset+1:offset+numFiles    %{'File','View','Move','Group','X','Y','Path'}
    [p,f,e] = fileparts(gd.DataSets(findex).filename);
    contents{findex, 1} = [f,e];
    contents{findex, 2} = false;
    contents{findex, 3} = false;
    contents{findex, 4} = gd.DataSets(findex).Map.YWorldLimits(1) - 0.5;
    contents{findex, 5} = gd.DataSets(findex).Map.XWorldLimits(1) - 0.5;
    contents{findex, 6} = p;
end
contents{1, 2} = true; % view first file
set(gd.Control.selection, 'Data', contents, 'Enable', 'on');
set(gd.Control.save, 'Enable', 'on');

guidata(hObject, gd);
plotDataAxes(gd);

% function gd = AdjustColorMaps(gd, indices)
% for findex = indices
%     gd.Experiment(findex).current = gd.Experiment(findex).Image;
%     % Adjust colormap
%     for cindex = 2
%         temp = gd.Experiment(findex).current(:,:,cindex);
%         temp = temp./gd.Experiment(findex).clim(cindex).current(2);
%         temp(temp>1) = 1;
%         temp(temp<gd.Experiment(findex).clim(cindex).current(1)/gd.Experiment(findex).clim(cindex).current(2)) = gd.Experiment(findex).clim(cindex).current(1)/gd.Experiment(findex).clim(cindex).current(2);
%         gd.Experiment(findex).current(:,:,cindex) = temp;
%     end
% end

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
        ImageA = gd.DataSets(viewSelection(fixed(1))).display;
    else
        ImageA = 1 - gd.DataSets(viewSelection(fixed(1))).display;
    end
    MapA = gd.DataSets(viewSelection(fixed(1))).Map;
    for index = 2:numel(fixed)
        if ~invert
            ImageA2 = gd.DataSets(viewSelection(fixed(index))).display;
        else
            ImageA2 = 1 - gd.DataSets(viewSelection(fixed(index))).display;
        end
        [ImageA, MapA] = imfuse(...
            ImageA,...
            MapA,...
            ImageA2,...
            gd.DataSets(viewSelection(fixed(index))).Map,...
            'blend',...
            'Scaling', 'Independent');
    end
else
    ImageA = [];
end
if ~isempty(moving)
    if ~invert
        ImageB = gd.DataSets(viewSelection(moving(1))).display;
    else
        ImageB = 1 - gd.DataSets(viewSelection(moving(1))).display;
    end
    MapB = gd.DataSets(viewSelection(moving(1))).Map;
    for index = 2:numel(moving)
        if ~invert
            ImageB2 = gd.DataSets(viewSelection(moving(index))).display;
        else
            ImageB2 = 1 - gd.gd.DataSets(viewSelection(moving(index))).display;
        end
        [ImageB, MapB] = imfuse(...
            ImageB,...
            MapB,...
            ImageB2,...
            gd.DataSets(viewSelection(moving(index))).Map,...
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

if eventdata.Indices(2) == 2 || eventdata.Indices(2) == 3       % display selected images
    plotDataAxes(gd);
    
elseif eventdata.Indices(2) == 4 || eventdata.Indices(2) == 5   % move image
    % Check input is numeric and valid
    if ~isnumeric(eventdata.NewData) || isnan(eventdata.NewData)
        hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
        return
    end
    % Record shift
    if eventdata.Indices(2) == 4    % up/down
        gd.DataSets(eventdata.Indices(1)).Map.YWorldLimits = [.5,gd.DataSets(eventdata.Indices(1)).Map.ImageSize(1)+.5] + eventdata.NewData;
    else                            % left/right
        gd.DataSets(eventdata.Indices(1)).Map.XWorldLimits = [.5,gd.DataSets(eventdata.Indices(1)).Map.ImageSize(2)+.5] + eventdata.NewData;
    end
    guidata(hObject, gd);   % save guidata
    plotDataAxes(gd);       % display shift
    
elseif eventdata.Indices(2) == 1 || eventdata.Indices(2) == 6   % adjust filename
    contents = get(hObject, 'Data'); 
    gd.DataSets(eventdata.Indices(1)).filename = fullfile(contents{eventdata.Indices(1), 6}, contents{eventdata.Indices(1), 1});
    guidata(hObject, gd);

elseif eventdata.Indices(2) == 7                    % remove selected file
    gd.DataSets(eventdata.Indices(1)) = [];         % remove from struct
    contents = get(hObject, 'Data');                % get table contents
    contents(eventdata.Indices(1),:) = [];          % remove from table
    set(hObject, 'Data', contents);                 % save table contents
    guidata(hObject, gd);
    plotDataAxes(gd);
end


function SaveData(hObject,eventdata,gd)
answer = questdlg('Are you sure? Any past Maps will be overwritten');
if ~strcmp(answer, 'Yes')
    return
end
set(gcf, 'pointer', 'watch')
drawnow
for findex = 1:numel(gd.DataSets);
    if exist(gd.DataSets(findex).filename, 'file')
        Map = gd.DataSets(findex).Map;
        save(gd.DataSets(findex).filename, 'Map', '-mat', '-append');
        fprintf('Map saved to: %s\n', gd.DataSets(findex).filename)
    else 
        
    end
end
set(gcf, 'pointer', 'arrow')

% function CloseReqFcn(src, ~)
% pause(.5);
% delete(gcf);
