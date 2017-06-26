function [saveFile, CLim] = vidDFoF(saveFile, Images, Maps, varargin)

StimIndex = [2 inf]; % stimuli indices to save to file
StimFrameIndex = [23, 47]; % frame indices of stimulus period for each stimulus
ControlIndex = 1; % false if no control trial
StimID = []; % ID #'s to display (including control trial)

% Display settings
type = 'AvgTrialdFoF'; % loading only
CMapType = 'parula';
showStimID = false;
showStimMarker = false;
showColorBar = true;
frameRate = 15.45*2;
mergetype = 'pretty'; % 'quick' or 'pretty'
Crop = false;
% Crop = [32.51, 0, 729.98, 512];
% Crop = [37.51,.51,722,510];
borderLims = [3,3,34,42];     % number of pixels to ignore from edges when computing ROI means, inclusive (top, bottom, left, right)
CLim = [];
filt = false;
% filt = fspecial('gaussian',5,1);
pixelSize = [1,1];
outputSize = [1080,1920/2];

directory = cd;

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'showStimID'
                showStimID = true;
                index = index + 1;
            case 'showStimMarker'
                showStimMarker = true;
                index = index + 1;
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
            case 'filter'
                filt = varargin{index+1};
                index = index + 2;
            case 'pixelSize'
                pixelSize = varargin{index+1};
                index = index + 2;
            case 'type'
                type = varargin{index+1};
                index = index + 2;
            case 'CMapType'
                CMapType = varargin{index+1};
                index = index + 2;
            case 'StimFrameIndex'
                StimFrameIndex = varargin{index+1};
                index = index + 2;
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case 'StimID'
                StimID = varargin{index+1};
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
    [Images,p] = uigetfile({'*.exp;*.align'}, 'Select image files for first dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(Images)
        return
    elseif iscellstr(Images)
        Images = fullfile(p, Images);
    else
        Images = {fullfile(p, Images)};
    end
    directory = p;
elseif ischar(Images)
    Images = {Images};
end

if ~exist('Maps', 'var') || isempty(Maps)
    [Maps,p] = uigetfile({'*.exp;*.align'}, 'Select map files for first dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(Maps)
        return
    elseif iscellstr(Maps)
        Maps = fullfile(p, Maps);
    else
        Maps = {fullfile(p, Maps)};
    end
elseif ischar(Maps)
    Maps = {Maps};
end


%% Load in images from the datasets
if iscellstr(Images)
    ImageFiles = Images;
    Images = cell(numel(ImageFiles),1);
    for findex = 1:numel(ImageFiles)
        temp = load(ImageFiles{findex}, type, '-mat');
        Images{findex} = temp.(type);
    end
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

% Set border to zero
if any(borderLims)
    for findex = 1:numel(Images)
        for sindex = 1:numel(Images{findex})
            Images{findex}{sindex}(1:borderLims(1),:,:,:,:) = 0;            % top
            Images{findex}{sindex}(end-borderLims(2)+1:end,:,:,:,:) = 0;    % bottom
            Images{findex}{sindex}(:,1:borderLims(3),:,:,:) = 0;            % left
            Images{findex}{sindex}(:,end-borderLims(4)+1:end,:,:,:) = 0;    % right
        end
    end
end

% Create map for pretty merge
if strcmp(mergetype, 'pretty')
    [Dim, Map, indMap] = mapFoVs(Maps, 'type', 'blend');
else
    indMap = [];
    Dim = [];
    Map = [];
end


%% Determine stimuli to save
if StimIndex(end) == inf
    StimIndex = [StimIndex(1:end-1), StimIndex(end-1)+1:numel(Images{1})];
end
numStims = numel(StimIndex);
if isempty(StimID)
    StimID = 1:numStims;
end

% Determine stimulus periods
if size(StimFrameIndex,1) == 1 && numel(StimIndex) > 1
    StimFrameIndex = repmat(StimFrameIndex,numStims,1);
end
   

%% Determine color info
if isempty(CLim)
    
    % Create A frame
    temp = cell(numFiles, 1);
    for ind = 1:numFiles
        temp{ind} = Images{ind}{ceil(end/2)}(:,:,StimFrameIndex(round(end/2),2));
    end
    Image = createImage(temp, Maps, 'speed', mergetype, 'filter', filt, 'Map', indMap, Dim, Map);

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
    case 'parula'
        cmap = parula(128);
end


%% Save each stimulus average to video

% Open video
fprintf('Writing video: %s...', saveFile);
vidObj = VideoWriter(saveFile,'Motion JPEG AVI');
set(vidObj, 'FrameRate', frameRate);
open(vidObj);

for index = 1:numStims
    sindex = StimIndex(index);
    
    for findex = 1:size(Images{1}{sindex},3)
        
        % Create frame
        temp = cell(numFiles, 1);
        for ind = 1:numFiles
            temp{ind} = Images{ind}{sindex}(:,:,findex);
        end
        Image = createImage(temp, Maps, 'speed', mergetype, 'filter', filt, 'Map', indMap, Dim, Map);
        
        % Scale for aspect ratio and desired output size
        sz = size(Image).*pixelSize; % determine how to scale to fix aspect ratio
        if ~isempty(outputSize)
            sz = sz*min(outputSize./sz); % determine conversion to reach desired size
        end
        Image = imresize(Image, sz); % scale video
        
        if ~isempty(outputSize) % add blanks to short edge to make output desired size
            sz = size(Image);
            dimind = find(sz-outputSize);
            num = outputSize(dimind)-sz(dimind);
            if dimind==1
                Image = cat(dimind,zeros(floor(num/2),sz(2)),Image,zeros(ceil(num/2),sz(2)));
            else
                Image = cat(dimind,zeros(sz(1),floor(num/2)),Image,zeros(sz(1),ceil(num/2)));
            end
            % Image = cat(2,Image,zeros(outputSize)); % SPECIFIC TO NEWS BLURB
        end
        
        % Display and save image to file
        if ~showStimID && ~showStimMarker && ~showColorBar % don't need to display image as no overlays
            Image = (Image-CLim(1))/range(CLim); % scale by CLim
            Image = round(Image*(size(cmap,1)-1)+1); % convert to colormap indexing
            Image = ind2rgb(Image, cmap);        % convert to RGB image
            writeVideo(vidObj, Image);           % save to file
            
        else % need to display image before saving
            
            if sindex == StimIndex(1) && findex == 1
                % Create figure
                hF = figure('Units', 'Pixels', 'Position', [50, 50, 1450, 950], 'Color', 'w');
                hA = axes('Parent', hF);
                
                % Determine frame dimensions
                [H, W, ~] = size(Image);
                StimH = round(H*pixelSize(1)/20);
                StimW = round(W*pixelSize(2)/20);
            end
            
            % Display Image
            imagesc(Image, CLim);
            axis off; hold on;
            colormap(cmap);
            
            % Place Stimulus mark
            if showStimMarker && sindex ~= ControlIndex && findex>=StimFrameIndex(index,1) && findex<=StimFrameIndex(index,2)
                patch([W-StimW*2; W-StimW; W-StimW; W-StimW*2],...
                    [H-StimH*2; H-StimH*2; H-StimH; H-StimH],...
                    'magenta','EdgeColor','magenta');
            end
            
            % Place ID number
            if showStimID
                if sindex == ControlIndex
                    text(StimH, StimW, 'control', 'FontSize', 20, 'Color', 'm', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                else
                    text(StimH, StimW, sprintf('%d', StimID(index)), 'FontSize', 20, 'Color', 'm', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
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
            if showColorBar
                frame = getframe(hF);
            else
                frame = getframe(hA);
            end
            writeVideo(vidObj, frame.cdata);
            
            % Close the figure
            if sindex == StimIndex(end) && findex == size(Images{1}{sindex},3)
                close(hF);
            end
            
        end % display
        
    end % frames

end % stims

close(vidObj);

fprintf('\tfinished\n');