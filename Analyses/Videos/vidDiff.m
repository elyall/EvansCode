function [saveFile, CLim] = vidDiff(saveFile, ImagesA, MapsA, ImagesB, MapsB, varargin)

StimIndex = [2 inf]; % stimuli indices to save to file
StimFrameIndex = [23, 47]; % frame indices of stimulus period for each stimulus
ControlIndex = 1; % false if no control trial
StimID = 0:8; % ID #'s to display (including control trial)

% Display settings
type = 'AvgTrialdFoF';
CMapType = 'HiLo';
showStimID = true;
showStimMarker = true;
showColorBar = true;
frameRate = 15.45*3;
mergetype = 'pretty'; % 'quick' or 'pretty'
% Crop = false;
Crop = [32.51, 0, 729.98, 512];
CLim = [];
filt = false;
% filt = fspecial('gaussian',5,1);

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
            case 'type'
                type = varargin{index+1};
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

if ~exist('ImagesA', 'var') || isempty(ImagesA)
    [ImagesA,p] = uigetfile({'*.exp;*.align'}, 'Select image files for first dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(ImagesA)
        return
    end
    ImagesA = fullfile(p, ImagesA);
    directory = p;
end
if ischar(ImagesA)
    ImagesA = {ImagesA};
end

if ~exist('MapsA', 'var') || isempty(MapsA)
    [MapsA,p] = uigetfile({'*.exp;*.align'}, 'Select map files for first dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(MapsA)
        return
    end
    MapsA = fullfile(p, MapsA);
end
if ischar(MapsA)
    MapsA = {MapsA};
end

if ~exist('ImagesB', 'var') || isempty(ImagesB)
    [ImagesB,p] = uigetfile({'*.exp;*.align'}, 'Select image files for second dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(ImagesB)
        return
    end
    ImagesB = fullfile(p, ImagesB);
    directory = p;
end
if ischar(ImagesB)
    ImagesB = {ImagesB};
end

if ~exist('MapsB', 'var') || isempty(MapsB)
    [MapsB,p] = uigetfile({'*.exp;*.align'}, 'Select map files for second dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(MapsB)
        return
    end
    MapsB = fullfile(p, MapsB);
end
if ischar(MapsB)
    MapsB = {MapsB};
end


%% Load in images from the two datasets
if iscellstr(ImagesA)
    ImageFiles = ImagesA;
    ImagesA = cell(numel(ImageFiles),1);
    for findex = 1:numel(ImageFiles)
        temp = load(ImageFiles{findex}, type, '-mat');
        ImagesA{findex} = temp.(type);
    end
end
numA = numel(ImagesA);

if iscellstr(ImagesB)
    ImageFiles = ImagesB;
    ImagesB = cell(numel(ImageFiles),1);
    for findex = 1:numel(ImageFiles)
        temp = load(ImageFiles{findex}, type, '-mat');
        ImagesB{findex} = temp.(type);
    end
end
numB = numel(ImagesB);


%% Load in maps from the two datasets
if iscellstr(MapsA)
    MapFiles = MapsA;
    MapsA = imref2d();
    for findex = 1:numel(MapFiles)
        load(MapFiles{findex}, 'Map', '-mat');
        MapsA(findex) = Map;
    end
end

if iscellstr(MapsB)
    MapFiles = MapsB;
    MapsB = imref2d();
    for findex = 1:numel(MapFiles)
        load(MapFiles{findex}, 'Map', '-mat');
        MapsB(findex) = Map;
    end
end

% Crop images
if ~islogical(Crop) || Crop ~= false
    [ImagesA, MapsA] = crop(ImagesA, Crop, MapsA);
    [ImagesB, MapsB] = crop(ImagesB, Crop, MapsB);
end

% Create map for pretty merge
if strcmp(mergetype, 'pretty')
    [DimA, MapA, indMapA] = mapFoVs(MapsA, 'type', 'blend');
    [DimB, MapB, indMapB] = mapFoVs(MapsB, 'type', 'blend');
else
    indMapA = [];
    MapA = [];
    DimA = [];
    indMapB = [];
    MapB = [];
    DimB = [];
end

% Initialize Map for Generating Images
[~, outMap] = mapFoVs([MapsA, MapsB]);
% outMap = mapOverlap(MapsA, MapsB, 'crop', Crop); % line above should do same thing


%% Determine stimuli to save
numStims = numel(ImagesA{1});
if StimIndex(end) == inf
    StimIndex = [StimIndex(1:end-1), StimIndex(end-1)+1:numStims];
end

% Determine stimulus periods
if size(StimFrameIndex,1) == 1 && numel(StimIndex) > 1
    StimFrameIndex = repmat(StimFrameIndex,numStims,1);
end


%% Determine color info
if isempty(CLim)
    
    % Create A frame
    temp = cell(numA, 1);
    for index = 1:numA
        temp{index} = ImagesA{index}{5}(:,:,1,StimFrameIndex(5,2));
    end
    ImageA = createImage(temp, MapsA, 'speed', mergetype, 'OutputView', outMap, 'filter', filt, 'Map', indMapA, DimA, MapA);
    
    % Create B frame
    temp = cell(numB, 1);
    for index = 1:numB
        temp{index} = ImagesB{index}{5}(:,:,1,StimFrameIndex(5,2));
    end
    ImageB = createImage(temp, MapsB, 'speed', mergetype, 'OutputView', outMap, 'filter', filt, 'Map', indMapB, DimB, MapB);
    
    % Compute difference
    Image = (ImageB - ImageA);%./(ImageB + ImageA);
        
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

% % stimmarker color
% if showStimMarker
%     cmap = cat(2, cmap, [1,0,1]);
%     showStimMarker = size(cmap,1);
% end


%% Save each stimulus average to video

% Open video
fprintf('Writing video: %s...', saveFile);
vidObj = VideoWriter(saveFile,'Motion JPEG AVI');
set(vidObj, 'FrameRate', frameRate);
open(vidObj);

% Create figure
hF = figure('Units', 'Pixels', 'Position', [50, 50, 1400, 900], 'Color', 'w');
hA = axes('Parent', hF);

for sindex = StimIndex
    
    for findex = 1:size(ImagesA{1}{sindex},4)
        
        % Create A frame
        temp = cell(numA, 1);
        for index = 1:numA
            temp{index} = ImagesA{index}{sindex}(:,:,1,findex);
        end
        ImageA = createImage(temp, MapsA, 'speed', mergetype, 'OutputView', outMap, 'filter', filt, 'Map', indMapA, DimA, MapA);
        
        % Create B frame
        temp = cell(numB, 1);
        for index = 1:numB
            temp{index} = ImagesB{index}{sindex}(:,:,1,findex);
        end
        ImageB = createImage(temp, MapsB, 'speed', mergetype, 'OutputView', outMap, 'filter', filt, 'Map', indMapB, DimB, MapB);
        
        % Compute difference
        Image = (ImageB - ImageA);%./(ImageB + ImageA);
        
        % Image = imfilter(Image, filt);
        
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
        if showStimID
            if sindex == ControlIndex
                text(StimH, StimW, 'control', 'FontSize', 20, 'Color', 'm', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
            else
                text(StimH, StimW, sprintf('%d', StimID(sindex)), 'FontSize', 20, 'Color', 'm', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
            end
        end
        
        % Display color bar
        if showColorBar
            cbH = colorbar;
            ylabel(cbH, 'Post (dF/F) - Pre (dF/F)');
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