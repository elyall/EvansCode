function [Filename, CLim, indMap, Dim, Map] = vidOverlay(Filename, Data, ColorIndex, varargin)
%vidOverlay Saves time-varying images of patch objects to a video file
%   FILENAME = vidOverlay() prompts user to select a FILENAME to save to,
%   and then select one or more ROI files to load ROIdata from.
%
%   FILENAME = vidOverlay(FILENAME, DATA, COLORINDEX) saves a video of DATA
%   to video file FILENAME. For accepted formats for DATA see OVERLAYROIS.
%   COLORINDEX is an array of size numROIs by numFrames assigning each ROI
%   a color for each frame by indexing into COLORS. COLORINDEX can also be
%   raw data that will be converted to indexing values.


% Default parameters that can be adjusted

% ROI properties
Maps = [];                  % cell array or array of imref2d objects specifying location of each input dataset
ROIindex = [];              % vector specifying ROIs to show
FileIndex = [];             % vector specifying corresponding file for each ROI
CLim = [];
Colors = 'b2r';             % Nx3 array of colors
roiType = 'roi';            % 'roi' or 'circle' specifying the type of ROI to plot
Radius = [];                % scalar greater than 0 specifying radius of circles
LineWidth = 1;              % scalar greater than 0 specifying width of outline
FaceAlpha = 1;              % scalar between 0 and 1 specifying opacity of ROI face
EdgeAlpha = 1;              % scalar between 0 and 1 specifying opacity of ROI outline
FaceBrightness = 1;         % scalar specifying brightness of ROI face
EdgeBrightness = 1;         % scalar specifying brightness of ROI edge

% Other overlays & colorbar
StimFrameIndex = [];        % frame indices or logical vector of length numFrames specifying which frames to display a stimulus marker on
Text = {};                  % cell array of text to overlay
TextIndex = [];             % array of length numFrames indexing which text to display on each frame
FontSize = 30;              % scalar specifying size of text
Color = [0,0,0];            % color of overlays
showColorbar = true;        % whether to display the colorbar
cbLabel = 'Fluorescence (A.U.)'; % string specifying ylabel on colorbar
cbFontSize = 20;            % scalar specifying size of colorbar text
flipColorbar = false;       % booleon specifying whether to flip colorbar up/down

% Display/video properties
Image = [];
CMap = [1,1,1];             % Nx3 colormap, or string specifying the colormap of the video
frameRate = 15.45*2;        % scalar specifying the frame rate of the output video
pixelSize = [1,1];          % two-element vector specifying the aspect ratio of the Y and X-axes respectively

% Placeholders
directory = cd; % default directory when prompting user to select a file

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Maps'
                Maps = varargin{index+1};
                index = index + 2;
            case 'ROIindex'
                ROIindex = varargin{index+1};
                index = index + 2;
            case 'FileIndex'
                FileIndex = varargin{index+1};
                index = index + 2;
            case 'CLim'
                CLim = varargin{index+1};
                index = index + 2;
            case 'Colors'
                Colors = varargin{index+1};
                index = index + 2;
            case 'roiType'
                roiType = varargin{index+1};
                index = index + 2;
            case 'Radius'
                Radius = varargin{index+1};
                index = index + 2;
            case 'LineWidth'
                LineWidth = varargin{index+1};
                index = index + 2;
            case 'FaceAlpha'
                FaceAlpha = varargin{index+1};
                index = index + 2;
            case 'EdgeAlpha'
                EdgeAlpha = varargin{index+1};
                index = index + 2;
            case 'FaceBrightness'
                FaceBrightness = varargin{index+1};
                index = index + 2;
            case 'EdgeBrightness'
                EdgeBrightness = varargin{index+1};
                index = index + 2;
            case 'StimFrameIndex'
                StimFrameIndex = varargin{index+1};
                index = index + 2;
            case 'Text'
                Text = varargin{index+1};
                index = index + 2;
            case 'TextIndex'
                TextIndex = varargin{index+1};
                index = index + 2;
            case 'FontSize'
                FontSize = varargin{index+1};
                index = index + 2;
            case 'Color'
                Color = varargin{index+1};
                index = index + 2;
            case {'Colorbar','colorbar','showColorbar'}
                showColorbar = varargin{index+1};
                index = index + 2;
            case 'cbLabel'
                cbLabel = varargin{index+1};
                index = index + 2;
            case 'cbFontSize'
                cbFontSize = varargin{index+1};
                index = index + 2;
            case 'flipColorbar'
                flipColorbar = true;
                index = index + 1;
            case 'Image'
                Image = varargin{index+1};
                index = index + 2;
            case {'CMap','Colormap','cmap','colormap'}
                CMap = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
                index = index + 2;
            case 'pixelSize'
                pixelSize = varargin{index+1};
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

if ~exist('Filename', 'var') || isempty(Filename)
    [Filename, p] = uiputfile({'*.avi'},'Save file as',directory);
    if isnumeric(Filename)
        return
    end
    Filename = fullfile(p, Filename);
end

if ~exist('Data','var') || isempty(Data)
    [Data, p] = uigetfile({'*.rois;*.mat'},'Choose ROI file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Data)
        return
    end
    Data = fullfile(p,Data);
end
if ischar(Data)
    Data = {Data};
end

if ~exist('ColorIndex','var') || isempty(ColorIndex)
    if iscell(Data) || isstruct(Data) || (iscell(Data) && isstruct(Data{1}))
        ColorIndex = gatherROIdata(Data,'dFoF',[],[],ROIindex,FileIndex);
        nR = size(ColorIndex,1);
        [nT,nF] = size(Data{1}.rois(1).dFoF);
        ColorIndex = reshape(ColorIndex,[nR,nT,nF]);
        ColorIndex = permute(ColorIndex,[1,3,2]);
        ColorIndex = reshape(ColorIndex,[nR,nT*nF]);
    else
        error('Specify color of ROIs input!');
    end
end
numFrames = size(ColorIndex,2);


%% Create background image

% Create image
if isempty(Image)
    [~, ~, indMap] = mapFoVs(Maps, 'type', 'index'); % generate indexing array
    Image = ones(size(indMap,1),size(indMap,2));
end

% Determine image colormap size
imgCLim = [min(Image(:)),max(Image(:))];
if isequal(round(Image),Image) && all(Image(:)>0) && range(imgCLim)+1<=256
    N = range(imgCLim)+1;
else
    N = 128;
end

% Create image colormap
if ischar(CMap)
    switch CMap
        case 'HiLo'
            CMap = HiLoColormap(imgCLim(1), imgCLim(2));
        case 'b2r'
            CMap = b2r(imgCLim(1),imgCLim(2));
        case 'gray'
            CMap = gray(N);
        case 'red'
            CMap = [linspace(0,1,N)', zeros(N,1), zeros(N,1)];
        case 'green'
            CMap = [zeros(N,1), linspace(0,1,N)', zeros(N,1)];
        case 'blue'
            CMap = [zeros(N,1), zeros(N,1), linspace(0,1,N)'];
        case 'white'
            CMap = repmat(linspace(0,1,N)',1,3);
        case 'parula'
            CMap = parula(N);
    end
end

% Convert image to index into colormap
if ~isequal(round(Image),Image) || any(Image(:)<1) || any(Image(:)>N)
    Image = scaleContinuousData(Image,[],'numSamples',N);
end


%% Load ROI colors

% Determine ROI colormap size
if isequal(round(ColorIndex),ColorIndex) && all(ColorIndex(:)>0) && range(CLim)+1<=256
    N = range(CLim)+1;
    if isempty(CLim)
        CLim = [min(ColorIndex(:)),max(ColorIndex(:))];
    end
else
    N = 128;
    if isempty(CLim)
        CLim = prctile(ColorIndex(:), [.01,99.99]);
    end
end

% Create ROI colormap
if ischar(Colors)
    switch Colors
        case 'HiLo'
            Colors = HiLoColormap(CLim(1),CLim(2));
            N = size(Colors,1);
        case 'b2r'
            Colors = b2r(CLim(1),CLim(2));
            N = size(Colors,1);
        case 'green'
            CMap = [zeros(N,1), linspace(0,1,N)', zeros(N,1)];
        case 'gray'
            Colors = gray(N);
        case 'parula'
            Colors = parula(N);
    end
end

% Convert data to index into colormap and determine tick marks
[ColorIndex, ~, ~, cbTick, cbTickLabel] = scaleContinuousData(ColorIndex,CLim,'numSamples',N);


%% Determine overlay info

% Add overlay color to colormap and determine overlay location
if any(StimFrameIndex) || any(TextIndex)
    [H, W, ~] = size(Image);
    Dist = min(round(H/10),round(W/10));
    Dist = Dist./pixelSize;
end

% Convert stim index to logical
if isnumeric(StimFrameIndex)
    temp = StimFrameIndex;
    StimFrameIndex = false(1,numFrames);
    StimFrameIndex(temp) = true;
end

% Conver numeric text to cellstr
if isnumeric(Text)
    if isrow(Text)
        Text = Text';
    end
    Text = cellstr(num2str(Text));
end

%% Display image and overlay ROIs

% Create figure
hF = figure('Units', 'Pixels', 'Position', [50, 50, 1450, 950], 'Color', 'w');
hA = axes('Parent', hF);

% Display Image and colorbar
image(Image);
colormap(CMap);
set(hA,'DataAspectRatio',[pixelSize 1]);
axis off; hold on;

% Overlay ROIs
patchHandles = overlayROIs(Data,...
    'axes',             hA,...
    'Maps',             Maps,...
    'ROIindex',         ROIindex,...
    'FileIndex',        FileIndex,...
    'roiType',          roiType,...
    'Radius',           Radius,...
    'Color',            Colors(ColorIndex(:,1),:),...
    'LineWidth',        LineWidth,...
    'FaceAlpha',        FaceAlpha,...
    'EdgeAlpha',        EdgeAlpha,...
    'FaceBrightness',   FaceBrightness,...
    'EdgeBrightness',   EdgeBrightness);
numROIs = numel(patchHandles);

% Display colorbar
if showColorbar
    
    % Create colorbar
    cbH = colorbar('westoutside','FontSize',cbFontSize/1.5);
    ylabel(cbH,cbLabel,'FontSize',cbFontSize);
    if flipColorbar
        set(cbH, 'Direction', 'reverse'); % flip colorbar
    end
    
    % Append colormap
    NewCMap = cat(1,CMap,Colors); % concatenate plot colors to colormap
    colormap(NewCMap);            % set new colormap
    
    % Determine limits of new section
    set(gca,'CLim',[0,size(NewCMap,1)]+.5);
    if ~flipColorbar
        NewYLim = [size(CMap,1),size(NewCMap,1)]+.5; % limit the colormap to this range
        cbTick = cbTick + size(CMap,1);
    else
        NewYLim = [0,size(Colors,1)]+.5;
    end
        
    % Center desired colormap and display labels
    set(cbH, 'Limits', NewYLim, 'FontSize', cbFontSize, 'Ticks', cbTick, 'YTickLabel', cbTickLabel);
    
end


%% Save data to video
fprintf('Writing video: %s...', Filename);

% Open video
vidObj = VideoWriter(Filename,'Motion JPEG AVI');
set(vidObj, 'FrameRate', frameRate);
open(vidObj);

% Save each frame to video
h = {};
axes(hA); % b2r bug
parfor_progress(numFrames);
for findex = 1:numFrames
    
    % Set ROI colors
    if findex~=1
        for rindex = 1:numROIs
            set(patchHandles(rindex),'FaceColor',Colors(ColorIndex(rindex,findex),:),...
                'EdgeColor',Colors(ColorIndex(rindex,findex),:));
        end
    end
    
    % Add stim marker
    if StimFrameIndex(findex)
        h{1} = patch([W-Dist(2)*2; W-Dist(2); W-Dist(2); W-Dist(2)*2],...
            [H-Dist(1)*2; H-Dist(1)*2; H-Dist(1); H-Dist(1)],Color,...
            'EdgeColor',Color);
    end
    
    % Add text
    if TextIndex(findex)
        h{2} = text(Dist(2)/2, Dist(1)/2, Text{TextIndex(findex)},'FontSize',FontSize,...
            'Color',Color,'HorizontalAlignment','left',...
            'VerticalAlignment','top');
    end
    
    % Write to video
    drawnow;
    frame = getframe(hF);
    Image = frame.cdata;
    writeVideo(vidObj, Image);
    
    % Remove frame-specific overlays
    if StimFrameIndex(findex)
        delete(h{1});
    end
    if TextIndex(findex)
        delete(h{2});
    end
    
    parfor_progress;
end

close(hF);
close(vidObj);
parfor_progress(0);

fprintf('\tfinished\n');

