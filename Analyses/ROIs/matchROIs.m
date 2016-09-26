function ROIindex = matchROIs(ROIs, Images, Maps, varargin)

ROIindex = [];

% default settings
gd.Internal.Settings.ROITextSize = 20;
gd.Internal.Settings.ROILineWidth = 2;


%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'ROIindex'
                ROIindex = varargin{index+1};
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
if ~exist('ROIs','var') || isempty(ROIs)
end
if ~iscell(ROIs)
    ROIs = {ROIs};
end
gd.Internal.numFiles = numel(ROIs);

if ~exist('Images','var') || isempty(Images)
end
if ~iscell(Images)
    Images = {Images};
end
if numel(Images) ~= gd.Internal.numFiles
    error('Number of Images needs to match number of ROI files input');
end

if ~exist('Maps','var') || isempty(Maps)
end
if numel(Maps) ~= gd.Internal.numFiles
    error('Number of Maps needs to match number of ROI files input');
end

%% Load in data
if iscellstr(ROIs)
    gd.Internal.ROIFiles = ROIs;
    for findex = 1:gd.Internal.numFiles
        load(ROIs{findex},'ROIdata','-mat');
        gd.ROIs{findex} = ROIdata; clear ROIdata;
    end
else
    gd.ROIs = ROIs;
end

if iscellstr(Images)
    gd.Internal.ImageFiles = Images;
    for findex = 1:gd.Internal.numFiles
        load(Images{findex},'m','-mat');
        gd.Images{findex} = m; clear m;
    end
else
    gd.Images = Images;
end

if iscellstr(Maps)
    gd.Internal.MapFiles = Maps;
    for findex = 1:gd.Internal.numFiles
        load(Maps{findex},'Map','-mat');
        gd.Maps(findex) = Map; clear Map;
    end
else
    gd.Maps = Maps;
end


%% Create ROIindex
gd.Internal.numROIs = cellfun(@(x) numel(x.rois), gd.ROIs);
if isempty(ROIindex)
    gd.ROIindex = [];
    for findex = 1:gd.Internal.numFiles
        temp = nan(gd.Internal.numROIs(findex),gd.Internal.numFiles);
        temp(:,findex) = 1:gd.Internal.numROIs(findex);
        gd.ROIindex = cat(1,gd.ROIindex,temp);
    end
else
    gd.ROIindex = ROIindex;
end
numROIs = size(gd.ROIindex,1);

% Generate colors
gd.Internal.colors.file = jet(gd.Internal.numFiles);
gd.Internal.colors.roi = jet(numROIs);
gd.Internal.colors.roi = gd.Internal.colors.roi(randperm(numROIs),:); % shuffle


%% Create & Populate Figure
% FIGURE
if nargout
    gd.fig = figure(...
        'NumberTitle',          'off',...
        'Name',                 'Match ROIs',...
        'ToolBar',              'figure',...
        'Units',                'pixels',...
        'Position',             [50, 50, 1400, 800],...
        'KeyPressFcn',          @(hObject,eventdata)KeyPressCallback(hObject,eventdata,guidata(hObject)),...
        'CloseRequestFcn',      'uiresume(gcbf)');
else
    gd.fig = figure(...
        'NumberTitle',          'off',...
        'Name',                 'Match ROIsw',...
        'ToolBar',              'figure',...
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
% axes 2
gd.Display.axes2 = axes(...
    'Parent',               gd.Display.panel,...
    'Units',                'normalized',...
    'Position',             [.5,0,.5,1],...
    'Visible',              'off');
axis off
% axes 1
gd.Display.axes1 = axes(...
    'Parent',               gd.Display.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,1,1]);
axis off
linkaxes([gd.Display.axes1,gd.Display.axes2]);

% CONTROLS
% panel
gd.Control.panel = uipanel(...
    'Title',                'Controls',...
    'Parent',               gd.fig,...
    'Units',                'Normalized',...
    'Position',             [.7,0,.3,1]);
% save button
gd.Control.save = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Save',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,.9,1,.1],...
    'Callback',             @(hObject,eventdata)SaveData(hObject,eventdata,guidata(hObject)),...
    'BackgroundColor',      [.7,1,.7],...
    'Enable',               'on');
% checkbox: merge or split
gd.Control.display = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Merge',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,.875,.25,.02],...
    'Value',                true,...
    'Callback',             @(hObject,eventdata)changePlot(hObject,eventdata,guidata(hObject)));
% checkbox: colors
gd.Control.color = uicontrol(...
    'Style',                'checkbox',...
    'String',               'Color Toggle',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [.25,.875,.25,.02],...
    'Value',                false,...
    'Callback',             @(hObject,eventdata)toggleColor(hObject,eventdata,guidata(hObject)));
% pushbutton: reset zoom
gd.Control.display = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Reset Zoom',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [.5,.875,.5,.025],...
    'Callback',             @(hObject,eventdata)resetZoom(hObject,eventdata,guidata(hObject)));
% table: cnames= name, view (logical), move (logical), side (popup: A or B)
gd.Control.data = uitable(...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,.775,1,.1],...
    'ColumnName',           {'View','View','Show All'},...
    'ColumnEditable',       [true,true,true],...
    'ColumnFormat',         {'logical','logical','logical'},...
    'ColumnWidth',          {100,100,100},...
    'Enable',               'on',...
    'CellEditCallback',     @(hObject,eventdata)DataSelection(hObject,eventdata,guidata(hObject)));
% button: merge
gd.Control.merge = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Merge',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,.675,.5,.1],...
    'Callback',             @(hObject,eventdata)MergeROIs(hObject,eventdata,guidata(hObject)),...
    'BackgroundColor',      [.7,1,.7],...
    'Enable',               'on');
% button: split
gd.Control.split = uicontrol(...
    'Style',                'pushbutton',...
    'String',               'Split',...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [.5,.675,.5,.1],...
    'Callback',             @(hObject,eventdata)SplitROIs(hObject,eventdata,guidata(hObject)),...
    'BackgroundColor',      [.7,1,.7],...
    'Enable',               'on');
% table: cnames= name, view (logical), move (logical), side (popup: A or B)
gd.Control.rois = uitable(...
    'Parent',               gd.Control.panel,...
    'Units',                'normalized',...
    'Position',             [0,0,1,.675],...
    'ColumnName',           repmat({'#','View'},1,gd.Internal.numFiles),...
    'ColumnEditable',       repmat([false,true],1,gd.Internal.numFiles),...
    'ColumnFormat',         repmat({'char','logical'},1,gd.Internal.numFiles),...
    'ColumnWidth',          repmat({30,40},1,gd.Internal.numFiles),...
    'Enable',               'on',...
    'CellEditCallback',     @(hObject,eventdata)ROISelection(hObject,eventdata,guidata(hObject))); %,...
    % 'CellSelectionCallback',@(src,evnt)set(src,'UserData',evnt.Indices));

% Create joint map
[gd.Map.dim, gd.Map.map, ~] = mapFoVs(gd.Maps,'type','index');

% Add files to table
gd.Control.data.Data = false(gd.Internal.numFiles,3); 
gd.Control.data.Data(1,1) = true;

% Add ROIs to table
temp = [num2cell(gd.ROIindex(:,1)),num2cell(false(size(gd.ROIindex,1),1))];
for findex = 2:gd.Internal.numFiles
    temp = [temp, [num2cell(gd.ROIindex(:,findex)),num2cell(false(size(gd.ROIindex,1),1))]];
end
gd.Control.rois.Data = temp;
gd.Internal.active = [];

gd.Internal.patch = []; % initialize handles matrix

guidata(gd.fig, gd);    % save guidata
plotDataAxes(gd);       % plot initial settings

% Return requested outputs
if nargout
    uiwait(gd.fig);
    gd = guidata(gd.fig);
    ROIindex = gd.ROIindex;
    delete(gd.fig)
end


%% CALLBACKS
function plotDataAxes(gd)

selection = gd.Control.data.Data(:,1:2);
offsets = gd.Map.dim;
currentZoom = [gd.Display.axes1.XLim; gd.Display.axes1.YLim];
if gd.Control.display.Value % Display one axis
    
    % Display only axis
    if any(selection(:,1))
        % offsets = mapFoVs(gd.Maps(selection(:,1)), 'Type', 'index');
        Image = createImage(gd.Images(selection(:,1)), gd.Maps(selection(:,1)),...
            'OutputView',   gd.Map.map);
        axes(gd.Display.axes1); hold off;
        h=imagesc(Image); axis off; hold on;
        set(h, 'ButtonDownFcn', @(hObject,eventdata)ROIClickCallback(hObject,eventdata,guidata(hObject)));
        if ~isequal(currentZoom,[0,1;0,1])
            gd.Display.axes1.XLim = currentZoom(1,:);
            gd.Display.axes1.YLim = currentZoom(2,:);
        end
    end
    
else % Display both axes
    
    % Display first axis
    if any(selection(:,1))
        Image = createImage(gd.Images(selection(:,1)), gd.Maps(selection(:,1)),...
            'OutputView',   gd.Map.map);
        axes(gd.Display.axes1); hold off;
        h=imagesc(Image); axis off; hold on;
        set(h, 'ButtonDownFcn', @(hObject,eventdata)ROIClickCallback(hObject,eventdata,guidata(hObject)));
        if ~isequal(currentZoom,[0,1;0,1])
            gd.Display.axes1.XLim = currentZoom(1,:);
            gd.Display.axes1.YLim = currentZoom(2,:);
        end
    end
    
    % Display second axis
    if any(selection(:,2))
        Image = createImage(gd.Images(selection(:,2)), gd.Maps(selection(:,2)),...
            'OutputView', gd.Map.map);
        axes(gd.Display.axes2); hold off;
        h=imagesc(Image); axis off; hold on;
        set(h, 'ButtonDownFcn', @(hObject,eventdata)ROIClickCallback(hObject,eventdata,guidata(hObject)));
        if ~isequal(currentZoom,[0,1;0,1])
            gd.Display.axes2.XLim = currentZoom(1,:);
            gd.Display.axes2.YLim = currentZoom(2,:);
        end
    end
end

% Display ROIs
index = cell2mat(gd.Control.rois.Data(:,2:2:end));
[rindex,findex] = find(index); % index in ROIindex
plotIndex = selection(findex,:);
if gd.Control.display.Value
    plotIndex(:,2) = false;
end
if ~isempty(rindex) && any(plotIndex(:))
    colors = jet(gd.Internal.numFiles);
    for r = 1:numel(rindex)
        if any(plotIndex(r,:))
            if plotIndex(r,1)
                axes(gd.Display.axes1);
                gd.Internal.patch(rindex(r),findex(r),1) = plotROI(rindex(r),findex(r),1,gd);
            end
            if plotIndex(r,2)
                axes(gd.Display.axes2);
                gd.Internal.patch(rindex(r),findex(r),2) = plotROI(rindex(r),findex(r),2,gd);
            end
        end
    end
    axes(gd.Display.axes1); hold off;
    if ~gd.Control.display.Value
        axes(gd.Display.axes2); hold off;
    end
end

guidata(gd.fig,gd); % update guidata

function h = plotROI(rindex,findex,axisindex,gd)
Rindex = gd.ROIindex(rindex,findex); % index in struct

% Translate vertices
vertices = bsxfun(@plus, gd.ROIs{findex}.rois(Rindex).vertices, gd.Map.dim(findex,1:2));

% Determine Linewidth
if ~isempty(gd.Internal.active) && ismember([rindex,findex], gd.Internal.active, 'rows')
    linewidth = gd.Internal.Settings.ROILineWidth*2;
else
    linewidth = gd.Internal.Settings.ROILineWidth;
end

if axisindex==1
    P = gd.Display.axes1;
else
    P = gd.Display.axes2;
end
if ~gd.Control.color.Value
    color = gd.Internal.colors.file(findex,:);
else
    color = gd.Internal.colors.roi(rindex,:);
end
% Plot ROI
h = patch(vertices(:,1),vertices(:,2),color,...
    'Parent',               P,...
    'FaceAlpha',            0,...
    'EdgeColor',            color,...
    'LineWidth',            linewidth,...
    'UserData',             [rindex,findex,Rindex],...
    'ButtonDownFcn',        @(hObject,eventdata)ROIClickCallback(hObject,eventdata,guidata(hObject))); % display ROI

function toggleColor(hObject,eventdata,gd)
[rindex,findex] = find(gd.Internal.patch);
if ~gd.Control.color.Value
    for f = 1:gd.Internal.numFiles
        set(gd.Internal.patch(rindex(findex==f),findex(findex==f)),'EdgeColor',gd.Internal.colors.file(f,:));
    end
else
    for r = 1:numel(rindex)
        
        set(gd.Internal.patch(rindex(r),findex(r)),'EdgeColor',gd.Internal.colors.roi(rindex(r),:));
    end
end

function changePlot(hObject,eventdata,gd)
if hObject.Value
    gd.Display.axes2.Visible = 'off';
    gd.Display.axes1.Position = [0,0,1,1]; axis off;
else
    gd.Display.axes1.Position = [0,0,.5,1]; axis off;
    gd.Display.axes2.Visible = 'on'; axis off;
end
plotDataAxes(gd);

function DataSelection(hObject,eventdata,gd)
if eventdata.Indices(2) == 3
    findex = eventdata.Indices(1);
    index = gd.ROIindex(:,findex)>0;
    if ~hObject.Data(eventdata.Indices(1),eventdata.Indices(2))
        [gd.Control.rois.Data{index,2*findex}] = deal(false); % display none
    else
        [gd.Control.rois.Data{index,2*findex}] = deal(true); % display all
    end
end
guidata(hObject,gd); % update guidata
plotDataAxes(gd);

function ROISelection(hObject,eventdata,gd)
index = [eventdata.Indices(1), eventdata.Indices(2), hObject.Data{eventdata.Indices(1),eventdata.Indices(2)-1}];
if ~isnan(index(3)) && ~isequal(index(3),0)   % current ROI exists
    if hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} % display ROI
        gd.Internal.patch(index(1),index(2)/2,1) = plotROI(index(1),index(2)/2,gd);
    else
        delete(gd.Internal.patch(index(1),index(2)/2)); % delete roi
        gd.Internal.patch(index(1),index(2)/2,1) = 0;
    end
else
    hObject.Data{eventdata.Indices(1),eventdata.Indices(2)} = false;
end
guidata(hObject,gd);

function ROIClickCallback(hObject,eventdata,gd)
if eventdata.Button == 1
    switch get(hObject, 'Type')
        case 'patch'
            [r,f] = find(gd.Internal.patch==hObject);
            set(gd.Internal.patch(r,f,1),'LineWidth',2*gd.Internal.Settings.ROILineWidth);
            gd.Internal.active = cat(1,gd.Internal.active,[r,f]);
        case 'image'
            for r = 1:size(gd.Internal.active,1)
                set(gd.Internal.patch(gd.Internal.active(r,1),gd.Internal.active(r,2),1),'LineWidth',gd.Internal.Settings.ROILineWidth);
            end
            gd.Internal.active = [];
    end
    guidata(hObject, gd);
end

function resetZoom(hObject,eventdata,gd)
gd.Display.axes1.XLim = gd.Map.map.XIntrinsicLimits;
gd.Display.axes1.YLim = gd.Map.map.YIntrinsicLimits;

function MergeROIs(hObject,eventdata,gd)
numActive = size(gd.Internal.active,1);
if numActive>1
    if numel(unique(gd.Internal.active(:,2)))<numActive
        error('Can''t link multiple ROIs from the same dataset!');
    end
    Rindex = min(gd.Internal.active(:,1));
    index = setdiff(gd.Internal.active(:,1),Rindex);
    [~,order] = sort(index,'descend');
    for r = order'
        attached = find(~isnan(gd.ROIindex(index(r),:)) & gd.ROIindex(index(r),:)~=0); % determine ROIs linked to current
        if any(~isnan(gd.ROIindex(Rindex,attached)) & gd.ROIindex(Rindex,attached)~=0) % make sure no ROIs are already linked in affected files
            error('Linked ROI already exists! Please unlink any already linked ROI');
        else
            
            % Update dictionary
            gd.ROIindex(Rindex,attached) = gd.ROIindex(index(r),attached);  % update registers
            gd.ROIindex(index(r),:) = [];                                   % remove old line from dictionary
            
            % Update display table
            for a = 1:numel(attached)
                gd.Control.rois.Data{Rindex,2*attached(a)-1} = gd.ROIindex(Rindex,attached(a)); % update table registers
                gd.Control.rois.Data{Rindex,2*attached(a)} = true;                              % update display properties
            end
            gd.Control.rois.Data(index(r),:) = []; % remove old line from table
            
            % Update patches and active list
            for a = 1:numel(attached)
                gd.Internal.patch(Rindex,attached(a),1) = gd.Internal.patch(index(r),attached(a),1); % transfer handle
                set(gd.Internal.patch(Rindex,attached(a),1),'UserData',[Rindex,attached(a),gd.ROIindex(Rindex,attached(a))]); % update userdata register
                if gd.Control.color.Value
                    set(gd.Internal.patch(Rindex,attached(a),1),'EdgeColor',gd.Internal.colors.roi(Rindex,:)); % update color
                end
            end
            gd.Internal.patch(index(r),:,:) = []; % remove old line from list
            gd.Internal.colors.roi(index(r),:) = [];  % remove old color

        end
    end
    
    % Update active list
    for a = gd.Internal.active(:,2)'
        set(gd.Internal.patch(Rindex,a,1),'LineWidth',gd.Internal.Settings.ROILineWidth); % remove from active
    end
    gd.Internal.active = []; % clear active list
    
    guidata(hObject, gd); % update guidata
end


function SplitROIs(hObject,eventdata,gd)
if size(gd.Internal.active,1) > 1
    error('Please only select one ROI for unlinking');
end
index = gd.Internal.active;
attached = find(~isnan(gd.ROIindex(index(1),:)) & gd.ROIindex(index(1),:)~=0); % determine ROIs linked to current
for r = 2:numel(attached)
    
    % Update dictionary
    gd.ROIindex(end+1,:) = nan;                                         % create new line
    gd.ROIindex(end,attached(r)) = gd.ROIindex(index(1),attached(r));   % update registers
    gd.ROIindex(index(1),attached(r)) = 0;                              %remove old register
    Rindex = size(gd.ROIindex,1);
    
    % Update display table
    [gd.Control.rois.Data{Rindex,2:2:end}] = deal(false);                       % create new line
    gd.Control.rois.Data{Rindex,2*attached(r)-1} = gd.ROIindex(end,attached(r));% update table register
    gd.Control.rois.Data{Rindex,2*attached(r)} = true;                          % update display
    gd.Control.rois.Data{index(1),2*attached(r)-1} = 0;                         % remove old registry from table
    gd.Control.rois.Data{index(1),2*attached(r)} = false;                       % update display

    % Update patch object
    gd.Internal.patch(Rindex,attached(r),1) = gd.Internal.patch(index(1),attached(r),1);  % transfer handle
    set(gd.Internal.patch(Rindex,attached(r),1),'UserData',[Rindex,attached(r),gd.ROIindex(end,attached(r))]); % update userdata register
    gd.Internal.colors.roi(Rindex,:) = rand(1,3); % generate new color
    set(gd.Internal.patch(Rindex,attached(r),1),'EdgeColor',gd.Internal.colors.roi(Rindex,:)); % update color
    gd.Internal.patch(index(1),attached(r),1) = 0; % remove old handle
    
    % Update active
    set(gd.Internal.patch(Rindex,attached(r),1),'LineWidth',gd.Internal.Settings.ROILineWidth);  % update linewidth
end
gd.Internal.active = []; % clear active list
set(gd.Internal.patch(index(1),attached(1),1),'LineWidth',gd.Internal.Settings.ROILineWidth);  % update linewidth
guidata(hObject, gd);


function SaveData(hObject,eventdata,gd)
[f,p] = uiputfile({'*.mat'},'Save ROIindex');
ROIindex = gd.ROIindex;
save(fullfile(p,f),'ROIindex','-mat');

