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
% gd.numFrames = size(Images, 4);
gd.numFrames = size(Images, 5);

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
gd.axes = axes();

% Create slider
minorstep = 1/(gd.numFrames-1);
gd.slider = uicontrol(...
    'Style',                'slider',...
    'Parent',               gd.fig,...
    'Units',                'normalized',...
    'Position',             [0,.95,1,.05],...
    'Min',                  1,...
    'Max',                  gd.numFrames,...
    'Value',                1,...
    'SliderStep',           [minorstep,max(2*minorstep,.1)],...
    'Callback',             @(hObject,eventdata)plotmainaxes(round(get(hObject, 'Value')), eventdata, guidata(hObject)));

guidata(gd.fig, gd);

% Plot first image
plotmainaxes(1, [], gd);

% If output: wait for figure to be closed
% waitfor(gd.fig);

function plotmainaxes(frameIndex, ~, gd)
axes(gd.axes)
% imagesc(gd.Images(:,:,1,frameIndex), gd.CLim)
imagesc(gd.Images(:,:,1,1,frameIndex), gd.CLim);
colormap(gd.CMap)
xlabel(sprintf('Frame %d', frameIndex));

