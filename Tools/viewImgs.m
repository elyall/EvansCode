function gd = viewImgs(Data,varargin)

CLim = [];
CMap = [];
FrameIndex = [1 inf]; % loading only
frameRate = 10;
summary = [];
pixels = [];
WhiskIDs = [];
points = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'CLim'
                CLim = varargin{index+1};
                index = index + 2;
            case 'CMap'
                CMap = varargin{index+1};
                index = index + 2;
            case {'Frames','frames'}
                FrameIndex = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
                index = index + 2;
            case 'summary'
                summary = varargin{index+1};
                index = index + 2;
            case 'pixels'
                pixels = varargin{index+1};
                index = index + 2;
            case 'WhiskIDs'
                WhiskIDs = varargin{index+1};
                index = index + 2;
            case 'points'
                points = varargin{index+1};
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

% Select and load in data
if ~exist('Data', 'var')
    [f,p] = uigetfile({'*.tif;*.sbx'}, 'Select image file');
    if isnumeric(f)
        return
    end
    Data = fullfile(p,f);
end
if ischar(Data)
    [~,~,ext] = fileparts(Data);
    switch ext
        case '.tif'
            gd.data.data = readTiff(Data, 'Frames', FrameIndex);
        case '.avi'
            gd.data.data = readAVI(Data, 'Frames', FrameIndex);
        case '.sbx'
            [gd.data.data, gd.data.loadObj] = load2P(Data, 'Type', 'Direct', 'Frames', FrameIndex, 'Double');
    end
else
    gd.data.data = Data;
end
gd.data.class = class(gd.data.data);
gd.data.shape = size(gd.data.data);

% Load in whisk variables
gd.data.whisk = false;
if ~isempty(summary)
    gd.data.whisk = true;
    if ischar(summary)
        [gd.data.summary,gd.data.pixels] = loadWhisk(summary,'CatOut');
    else
        gd.data.summary = summary;
        gd.data.pixels = pixels;
    end
    if isempty(WhiskIDs)
        [temp,~,gd.data.ids] = unique(gd.data.summary.label);
    else
        [temp,~,gd.data.ids] = unique(WhiskIDs);
    end
    gd.data.colors = lines(numel(temp));
end

% Load in DLC variables
gd.data.DLC = false;
if ~isempty(points)
    gd.data.DLC = true;
    if ischar(points)
        gd.data.xy = loadDLC(points);
    else
        gd.data.xy = points;
    end
    gd.data.DLCcolors = jet(5);
end

% Determine color limits
if isempty(CLim)
%     gd.internal.CLim = prctile(Images(linspace(1,numel(Images),prod(gd.data.shape(1:2))*100)), [1,99]);
    gd.internal.CLim = double([min(gd.data.data(:)), max(gd.data.data(:))]);
else
    gd.internal.CLim = CLim;
end

% Save colormap input
if isempty(CMap)
    gd.internal.CMap = [];
else
    gd.internal.CMap = CMap;
end

%% Generate GUI

% Create figure
gd.gui.fig = figure(...
    'KeyPressFcn',      @(hObject,eventdata)KeyPressCallback(hObject,eventdata,guidata(hObject)));

% Create axis
gd.gui.axes = axes(...
    'Parent',           gd.gui.fig,...
    'Units',            'normalized',...
    'XLimMode',         'manual',...
    'YLimMode',         'manual');
if isequal(gd.data.class, 'logical')
    gd.gui.axes.Position = [.1, .1, .8, .7];
else
    gd.gui.axes.Position = [.1, .2, .8, .6];
end

% Create indexing sliders
maxvalue = gd.data.shape(4);
minorstep = 1/(maxvalue-1);
ToolTipString = {'frame','first','last'};
Values = [1,1,maxvalue];
for sindex = 1:3
    yloc = .8+(3-sindex)*.2/3;
    gd.gui.sliders(sindex) = uicontrol(...
        'Style',                'slider',...
        'Parent',               gd.gui.fig,...
        'Units',                'normalized',...
        'Position',             [.1,yloc,.9,.2/3],...
        'Min',                  1,...
        'Max',                  maxvalue,...
        'Value',                Values(sindex),...
        'SliderStep',           [minorstep,max(2*minorstep,.1)],...
        'ToolTipString',        ToolTipString{sindex},...
        'UserData',             sindex,...
        'Callback',             @(hObject,eventdata)updateindex(hObject, eventdata, guidata(hObject)),...
        'ButtonDownFcn',        @(hObject,eventdata)updateindex(hObject, eventdata, guidata(hObject)));
    gd.gui.sliderText(sindex) = uicontrol(...
        'Style',                'text',...
        'Parent',               gd.gui.fig,...
        'Units',                'normalized',...
        'Position',             [0,yloc,.1,.2/3],...
        'String',               sprintf('%d/%d',Values(sindex),maxvalue),...
        'HorizontalAlignment',  'right');
end

% Create video button
gd.gui.play = uicontrol(...
    'Style',            'togglebutton',...
    'Parent',           gd.gui.fig,...
    'String',           'Play',...
    'Interruptible',    'on',...
    'Units',            'normalized',...
    'Position',         [0,.05,.2,.05],...
    'Callback',         @(hObject,eventdata)Play(hObject,eventdata,guidata(hObject)));

% Create frame rate  input
gd.gui.frameRate = uicontrol(...
    'Style',            'edit',...
    'Parent',           gd.gui.fig,...
    'String',           frameRate,...
    'Units',            'normalized',...
    'Position',         [0,.1,.2,.05],...
    'Callback',         @(hObject,eventdata)FrameRate(hObject,eventdata,guidata(hObject)));

% Create video button
gd.gui.save = uicontrol(...
    'Style',            'togglebutton',...
    'Parent',           gd.gui.fig,...
    'String',           'Save',...
    'Interruptible',    'on',...
    'Units',            'normalized',...
    'Position',         [0,.15,.2,.05],...
    'Callback',         @(hObject,eventdata)Save(hObject,eventdata,guidata(hObject)));

% Create colormap dropdown
gd.gui.colormap = uicontrol(...
    'Style',            'popupmenu',...
    'Parent',           gd.gui.fig,...
    'String',           {'gray','parula'},...
    'Units',            'normalized',...
    'Position',         [0,0,.2,.05],...
    'Callback',         @(hObject,eventdata)plotmainaxes(guidata(hObject)));
if ~isempty(gd.internal.CMap)
    gd.gui.colormap.String{3} = 'input';
end

% Create colorlimit sliders
for index = 1:2
    minorstep = min(1/(gd.internal.CLim(2)-gd.internal.CLim(1)),.05);
    yloc = (index-1)*.1/2;
    gd.gui.clim(index) = uicontrol(...
        'Style',                'slider',...
        'Parent',               gd.gui.fig,...
        'Units',                'normalized',...
        'Position',             [.35,yloc,.45,.1/2],...
        'Min',                  gd.internal.CLim(1),...
        'Max',                  gd.internal.CLim(2),...
        'Value',                gd.internal.CLim(index),...
        'SliderStep',           [minorstep,max(2*minorstep,.1)],...
        'UserData',             index,...
        'Callback',             @(hObject,eventdata)updateCLim(hObject, eventdata, guidata(hObject)));
    gd.gui.climText(index) = uicontrol(...
        'Style',                'text',...
        'Parent',               gd.gui.fig,...
        'Units',                'normalized',...
        'Position',             [.2,yloc,.15,.1/2],...
        'String',               sprintf('%d', gd.internal.CLim(index)),...
        'HorizontalAlignment',  'right');
end

% Create whisker overlay toggle
if gd.data.whisk || gd.data.DLC
    gd.gui.whiskOverlay = uicontrol(...
        'Style',                'checkbox',...
        'Parent',               gd.gui.fig,...
        'Units',                'normalized',...
        'Position',             [.8,.1,.2,.05],...
        'Value',                1,...
        'String',               'overlay',...
        'Callback',             @(hObject,eventdata)plotmainaxes(guidata(hObject)));
end

% Create colorlimit individual toggle
gd.gui.climEach = uicontrol(...
    'Style',                'checkbox',...
    'Parent',               gd.gui.fig,...
    'Units',                'normalized',...
    'Position',             [.8,.05,.2,.05],...
    'Value',                0,...
    'String',               'CLim each',...
    'Callback',             @(hObject,eventdata)CLimEach(hObject, eventdata, guidata(hObject)));

% Create histeq toggle
gd.gui.histeq = uicontrol(...
    'Style',                'checkbox',...
    'Parent',               gd.gui.fig,...
    'Units',                'normalized',...
    'Position',             [.8,0,.2,.05],...
    'Value',                0,...
    'String',               'hist eq',...
    'Callback',             @(hObject,eventdata)plotmainaxes(guidata(hObject)));

guidata(gd.gui.fig,gd);

%% Initialization

plotmainaxes(gd); % plot first image

% Play video
% gd.gui.play.Value = true;
% Play(gd.gui.play,[],gd);

end

%% Key Press Callback
function KeyPressCallback(hObject, eventdata, gd)
switch eventdata.Key
    case 'space'
        gd.gui.play.Value = ~gd.gui.play.Value;
        Play(gd.gui.play,[],gd);
    case 'uparrow'
        gd.gui.frameRate.String = string(min(str2double(gd.gui.frameRate.String)+5,100));
    case 'downarrow'
        gd.gui.frameRate.String = string(max(str2double(gd.gui.frameRate.String)-5,0.001));
    case 'rightarrow'
        if gd.gui.sliders(1).Value == gd.data.shape(end)
            gd.gui.sliders(1).Value = gd.gui.sliders(2).Value;
        else
            gd.gui.sliders(1).Value = gd.gui.sliders(1).Value + 1;
        end
        updateindex(gd.gui.sliders(1),[],gd);
    case 'leftarrow'
        if gd.gui.sliders(1).Value == 1
            gd.gui.sliders(1).Value = gd.gui.sliders(3).Value;
        else
            gd.gui.sliders(1).Value = gd.gui.sliders(1).Value - 1;
        end
        updateindex(gd.gui.sliders(1),[],gd);
end
end

%% Adjust controls

function updateindex(hObject, ~, gd)
% Round values
for index = 1:3
    gd.gui.sliders(index).Value = round(gd.gui.sliders(index).Value);
end
% Ensure boundaries don't conflict
if gd.gui.sliders(2).Value>gd.gui.sliders(3).Value
    gd.gui.sliders(2).Value = gd.gui.sliders(3).Value;
end
if gd.gui.sliders(3).Value<gd.gui.sliders(2).Value
    gd.gui.sliders(3).Value = gd.gui.sliders(2).Value;
end
if hObject.UserData == 1
    % If selection is out of bounds, move to other side
    if gd.gui.sliders(1).Value<gd.gui.sliders(2).Value
        gd.gui.sliders(1).Value = gd.gui.sliders(3).Value;
    elseif gd.gui.sliders(1).Value>gd.gui.sliders(3).Value
        gd.gui.sliders(1).Value = gd.gui.sliders(2).Value;
    end
else
    % If selection is out of bounds, move inbounds
    if gd.gui.sliders(1).Value<gd.gui.sliders(2).Value
        gd.gui.sliders(1).Value = gd.gui.sliders(2).Value;
    elseif gd.gui.sliders(1).Value>gd.gui.sliders(3).Value
        gd.gui.sliders(1).Value = gd.gui.sliders(3).Value;
    end
end
% Update GUI
for index = 1:3
    gd.gui.sliderText(index).String = sprintf('%d/%d', gd.gui.sliders(index).Value, gd.data.shape(end)); % update text
end
plotmainaxes(gd) % plot new image
end

function Play(hObject,eventdata,gd)
if hObject.Value
    set(hObject,'String','Stop','BackgroundColor',[0,0,0],'ForegroundColor',[1,1,1]);
    while hObject.Value
        tic
        if gd.gui.sliders(1).Value == gd.data.shape(end)
            gd.gui.sliders(1).Value = gd.gui.sliders(2).Value;
        else
            gd.gui.sliders(1).Value = gd.gui.sliders(1).Value + 1;
        end
        updateindex(gd.gui.sliders(1),[],gd);
        pause(.007)
        while toc < 1/str2double(gd.gui.frameRate.String)
            pause(.007)
        end
    end
    set(hObject,'String','Play','BackgroundColor',[.94,.94,.94],'ForegroundColor',[0,0,0]);
else
    hObject.String = 'Stopping...';
end
end

function Save(hObject,eventdata,gd)
if hObject.Value
    set(hObject,'String','Saving...');
    [f,p] = uiputfile({'*.avi'},'Save video as');
    if isnumeric(f)
        set(hObject,'String','Saving...','BackgroundColor',[.94,.94,.94],'ForegroundColor',[0,0,0]);
        return
    end
    vidObj = VideoWriter(fullfile(p,f),'Motion JPEG AVI'); % initiate file
    set(vidObj, 'FrameRate', str2double(gd.gui.frameRate.String)); % set frame rate
    open(vidObj); % open file
    findex = gd.gui.sliders(2).Value-1;
    while hObject.Value && findex < gd.gui.sliders(3).Value
        findex = findex + 1;                  % increment index
        gd.gui.sliders(1).Value = findex;     % update slider
        updateindex(gd.gui.sliders(1),[],gd); % plot current frame
        frame = getframe(gd.gui.axes);        % get image
        writeVideo(vidObj,frame.cdata);       % write to video
    end
    close(vidObj);
    set(hObject,'String','Save');
    if hObject.Value
        set(hObject,'Value',0);
    end
end
end

function FrameRate(hObject,eventdata,gd)
Value = str2double(hObject.String);
if isnan(Value)
    hObject.String = 30;
elseif Value>100
    hObject.String = 100;
elseif Value<=0
    hObject.String = 0.001;
end
end

function updateCLim(hObject, ~, gd)
% update CLim
value = hObject.Value;
if hObject.UserData==1
    if value >= gd.internal.CLim(2)
        hObject.Value = hObject.Min;
        return
    end
elseif hObject.UserData==2
    if value <= gd.internal.CLim(1)
        hObject.Value = hObject.Max;
        return
    end
end
gd.internal.CLim(hObject.UserData) = value;  % set value
guidata(hObject, gd);               % save guidata
gd.gui.climText(hObject.UserData).String = sprintf('%d', gd.internal.CLim(hObject.UserData)); % update text
plotmainaxes(gd);                   % plot new image
end

function CLimEach(hObject, ~, gd)
if hObject.Value
    set([gd.gui.clim,gd.gui.climText], 'Enable', 'off'); % disable clim sliders
else
    set([gd.gui.clim,gd.gui.climText], 'Enable', 'on');  % enable clim sliders
end
plotmainaxes(gd); % plot new image
end


%% Display image

function plotmainaxes(gd)
axes(gd.gui.axes)
val = round(gd.gui.sliders(1).Value);
img = gd.data.data(:,:,1,val);
gd.gui.histeq.Enable = 'on';
if gd.gui.histeq.Value
    img = adapthisteq(img, 'NumTiles', [16 16], 'Distribution', 'Exponential');
end
switch gd.data.class
    case 'logical'
        imshow(img);
    otherwise
        if gd.gui.climEach.Value
            imagesc(img);
        else
            imagesc(img, gd.internal.CLim);
        end
end
axis equal off

% Set colormap
if gd.gui.colormap.Value==1
    colormap('gray')
elseif gd.gui.colormap.Value==2
    colormap('parula')
elseif gd.gui.colormap.Value==3
    colormap(gd.internal.CMap)
end

% Overlay whiskers
if gd.data.whisk && gd.gui.whiskOverlay.Value
    hold on;
    for n = find(gd.data.summary.time==val-1)'
%         plot(gd.data.summary.fol_x(n),gd.data.summary.fol_y(n),'bx');
        if ismember(gd.data.summary.label(n),gd.data.ids)
            plot(gd.data.pixels{n}(:,1),gd.data.pixels{n}(:,2),'-','Color',gd.data.colors(gd.data.ids(n),:)); % colored by ID
            text(double(gd.data.summary.fol_x(n)),double(gd.data.summary.fol_y(n)),sprintf('%d',gd.data.summary.label(n)),'Color',[1,1,1]);
%         else
%             plot(gd.data.pixels{n}(:,1),gd.data.pixels{n}(:,2),'r-'); % all same color
        end
    end
    hold off;
end

if gd.data.DLC && gd.gui.whiskOverlay.Value
    hold on;
    for n = 1:5
        plot(squeeze(gd.data.xy(val,1,3*(n-1)+1:3*n)),squeeze(gd.data.xy(val,2,3*(n-1)+1:3*n)),'x','Color',gd.data.DLCcolors(n,:));
    end
    hold off;
end

end

