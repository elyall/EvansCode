function displayStack(Images, CLim, CMap)


%% Parse input arguments

% Load images
if ischar(Images)
    ImageFiles = Images;
    load(ImageFiles, 'data');
    Images = data;
    [gd.Images, gd.loadObj] = load2P(ImageFiles, 'Type', 'Direct');
else
    gd.Images = Images;
end
gd.class = class(Images);

% Format images
while ndims(gd.Images) < 5
    gd.Images = permute(gd.Images, [1:ndims(gd.Images)-1, ndims(gd.Images)+1, ndims(gd.Images)]);
end
gd.dim = size(gd.Images);

% Determine color limits
if ~exist('CLim', 'var') || isempty(CLim)
%     gd.CLim = prctile(Images(linspace(1,numel(Images),512*796*100)), [1,99]);
    gd.CLim = [min(Images(:)), max(Images(:))];
else
    gd.CLim = CLim;
end

% Determine color map
if ~exist('CMap', 'var') || isempty(CMap)
    gd.CMap = 'gray';
%     gd.CMap = HiLoColormap(CLim(1), CLim(2));
else
    gd.CMap = CMap;
end

%% Code

% Create figure
gd.fig = figure;

% Create axes
gd.axes = axes('Parent', gd.fig, 'Units', 'normalized', 'Position', [0, 0, 1, .9]);

% Create sliders
sliderIndex = find(gd.dim(3:end)>1);
ndim = numel(sliderIndex);
for sindex = 1:ndim
    maxvalue = gd.dim(sliderIndex(sindex)+2);
    minorstep = 1/(maxvalue-1);
    yloc = .9+(ndim-sindex)*.1/ndim;
    gd.sliders(sindex) = uicontrol(...
        'Style',                'slider',...
        'Parent',               gd.fig,...
        'Units',                'normalized',...
        'Position',             [.1,yloc,.9,.1/ndim],...
        'Min',                  1,...
        'Max',                  maxvalue,...
        'Value',                1,...
        'SliderStep',           [minorstep,max(2*minorstep,.1)],...
        'ToolTipString',        sprintf('Dimension %d', sliderIndex(sindex)),...
        'UserData',             [sindex, sliderIndex(sindex), maxvalue],...
        'Callback',             @(hObject,eventdata)updateindex(hObject, eventdata, guidata(hObject)));
    gd.text(sindex) = uicontrol(...
        'Style',                'text',...
        'Parent',               gd.fig,...
        'Units',                'normalized',...
        'Position',             [0,yloc,.1,.1/ndim],...
        'String',               sprintf('1/%d', maxvalue),...
        'HorizontalAlignment',  'right');
end
gd.Position = ones(1,3);

guidata(gd.fig, gd);

% Plot first image
plotmainaxes(gd);

% If output: wait for figure to be closed
% waitfor(gd.fig);

function updateindex(hObject, ~, gd)
% update index
gd.Position(hObject.UserData(2)) = round(hObject.Value);
guidata(hObject, gd);

% update text
gd.text(hObject.UserData(1)).String = sprintf('%d/%d', gd.Position(hObject.UserData(2)), hObject.UserData(3));

% plot new image
plotmainaxes(gd)

function plotmainaxes(gd)
axes(gd.axes)
switch gd.class
    case 'logical'
        imshow(gd.Images(:,:, gd.Position(1), gd.Position(2), gd.Position(3)));
    otherwise
        imagesc(gd.Images(:,:, gd.Position(1), gd.Position(2), gd.Position(3)), gd.CLim);
end
colormap(gd.CMap)
axis off
% xlabel(sprintf('Frame %d', Index));
