function [saveFile, CLim] = vidRaw(saveFile, Images, Maps, varargin)

% Loading only (if cellstr or char input for Images)
Frames = [1 inf];
Channel = 1;

% Stimulus indicator (vector)
StimFrameIndex = false;

% Display settings
CMapType = 'green';
showColorBar = true;
frameRate = 15.45;
mergetype = 'pretty'; % 'quick' or 'pretty'
Crop = false;
% Crop = [32.51, 0, 729.98, 512];
CLim = [];
Afilt = false;
% Afilt = fspecial('gaussian',5,1);
Tfilt = false;
% Tfilt = fspecial('gaussian',5,1);
% Tfilt = Tfilt(round(5/2),:);

directory = cd;

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'showColorBar'
                showColorBar = true;
                index = index + 1;
            case 'Crop'
                Crop = varargin{index+1};
                index = index + 2;
            case 'CLim'
                CLim = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
                index = index + 2;
            case 'mergetype'
                mergetype = varargin{index+1};
                index = index + 2;
            case 'Afilter'
                Afilt = varargin{index+1};
                index = index + 2;
            case 'CMapType'
                CMapType = varargin{index+1};
                index = index + 2;
            case 'StimFrameIndex'
                StimFrameIndex = varargin{index+1};
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

if ~exist('saveFile', 'var') || isempty(saveFile)
    [saveFile, p] = uiputfile({'*.avi'},'Save file as',directory);
    if isnumeric(saveFile)
        return
    end
    saveFile = fullfile(p, saveFile);
end

if ~exist('Images', 'var') || isempty(Images)
    [Images,p] = uigetfile({'*.sbx;*.tiff'}, 'Select image files:', directory, 'MultiSelect', 'on');
    if isnumeric(Images)
        return
    end
    Images = fullfile(p,Images);
    directory = p;
elseif ischar(Images)
    Images = {Images};
end

if ~exist('Maps', 'var') || isempty(Maps)
    [Maps,p] = uigetfile({'*.exp;*.align'}, 'Select map files:', directory, 'MultiSelect', 'on');
    if isnumeric(Maps)
        return
    end
    Maps = fullfile(p,Maps);
elseif ischar(Maps)
    Maps = {Maps};
end


%% Load in images from the datasets
if iscellstr(Images)
    ImageFiles = Images;
    Images = cell(numel(ImageFiles),1);
    for findex = 1:numel(ImageFiles)
        Images{findex} = load2P(ImageFiles{findex},'Frames',Frames,'Channel',Channel);
    end
elseif isnumeric(Images)
    Images = {Images};
end
numFiles = numel(Images);


%% Load in maps from the datasets
if iscellstr(Maps)
    MapFiles = Maps;
    Maps = imref2d();
    for findex = 1:numel(MapFiles)
        load(MapFiles{findex}, 'Map', '-mat');
        Maps(findex) = Map;
    end
end

% Crop images
if ~islogical(Crop) || Crop ~= false
    [Images, Maps] = crop(Images, Crop, Maps);
end

% Create map for pretty merge
if strcmp(mergetype, 'pretty')
    [Dim, Map, indMap] = mapFoVs(Maps, 'type', 'blend');
else
    indMap = [];
    Dim = [];
    Map = [];
end
   

%% Determine color info
if isempty(CLim)
    
    % Determine Color Limits
    temp = cell(numFiles, 1);
    nFrames = cellfun(@(x) size(x,5), Images);
    for index = 1:numFiles
        temp{index} = Images{index}(:,:,1,round(nFrames(index)/2));
    end
    Image = createImage(temp, Maps, 'speed', mergetype, 'filter', Afilt, 'Map', indMap, Dim, Map);
    CLim = prctile(Image(:), [.01,99.99]);
    
end

% Determine colormap
switch CMapType
    case 'HiLo'
        cmap = HiLoColormap(CLim(1), CLim(2));
    case 'b2r'
        cmap = b2r(CLim(1),CLim(2));
    case 'gray'
        cmap = gray(128);
    case 'red'
        cmap = [linspace(0,1,128)', zeros(128,1), zeros(128,1)];
    case 'green'
        cmap = [zeros(128,1), linspace(0,1,128)', zeros(128,1)];
    case 'blue'
        cmap = [zeros(128,1), zeros(128,1), linspace(0,1,128)'];
end


%% Save each stimulus average to video

% Open video
fprintf('Writing video: %s...', saveFile);
vidObj = VideoWriter(saveFile,'Motion JPEG AVI');
set(vidObj, 'FrameRate', frameRate);
open(vidObj);

% Create figure
hF = figure('Units', 'Pixels', 'Position', [50, 50, 1450, 950], 'Color', 'w');
hA = axes('Parent', hF);

for sindex = StimIndex
    
    for findex = 1:size(Images{1}{sindex},4)
        
        % Create frame
        temp = cell(numFiles, 1);
        for index = 1:numFiles
            temp{index} = Images{index}{sindex}(:,:,1,findex);
        end
        Image = createImage(temp, Maps, 'speed', mergetype, 'filter', Afilt, 'Map', indMap, Dim, Map);
                
        % Display Image
        imagesc(Image, CLim);
        axis off;
        colormap(cmap);
        
        % Determine frame dimensions
        if sindex == StimIndex(1) && findex == 1
            [H, W, ~] = size(Image);
            StimH = round(H/20);
            StimW = round(W/20);
            hold on;
        end
        
        % Place Stimulus mark
        if showStimMarker && sindex ~= ControlIndex && findex>=StimFrameIndex(sindex,1) && findex<=StimFrameIndex(sindex,2)
            patch([W-StimW*2; W-StimW; W-StimW; W-StimW*2],...
                [H-StimH*2; H-StimH*2; H-StimH; H-StimH],...
                'magenta','EdgeColor','magenta');
        end
        
        % Place ID number
        if StimIndex
            if sindex == ControlIndex
                text(StimH, StimW, 'control', 'FontSize', 20, 'Color', 'm', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
            else
                text(StimH, StimW, sprintf('%d', StimID(sindex)), 'FontSize', 20, 'Color', 'm', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
            end
        end
        
        % Display color bar
        if showColorBar
            cbH = colorbar;
            ylabel(cbH, 'dF/F');
        end
        
        % Write to video
        drawnow
        pause(0.01);
        frame = getframe(hF);
        writeVideo(vidObj, frame.cdata);
        
    end

end

close(vidObj);
close(hF);

fprintf('\tfinished\n');