function [Filename, CLim] = vidDFoF(Filename, Images, Maps, varargin)
%vidDFoF Saves images to a video file
%   SAVEFile = vidDFoF() prompts user to select a FILENAME to save to and
%   select a .sbx or .tif file to load data from.

% TO DO: use insertText to display StimIDs without displaying to figure.


% Default parameters that can be adjusted

ControlIndex = 1;           % index in cell array, false if no control trial
StimFrameIndex = [];        % frame indices of stimulus period for each stimulus
StimID = {};                % text to display for each stimuls
VarToLoad = 'AvgTrialdFoF'; % variable name to load from IMAGES if IMAGES is a char or cell array of strings
CMap = 'parula';            % colormap or string specifying which colormap to use
CLim = [];                  % 1x2 vector specifying the color limits ([] uses max and min)
mergetype = 'pretty';       % 'quick' or 'pretty'
filt = false;               % 2D filter to apply to each frame
% filt = fspecial('gaussian',5,1);               % 2D filter to apply to each frame
showColorBar = false;       % specifies whether to display the colorbar
frameRate = 15.45*2;        % specifies the frame rate of the output video
Crop = false;               % number of pixels to remove from edges (see: crop)
borderLims = [];            % number of pixels along edges to set to zero, inclusive (top, bottom, left, right)
pixelSize = [1,1];          % 1x2 vector specifying the aspect ratio (um/pixel)
outputSize = [];            % 1x2 vector specifying the desired Height and Width of the output video


% Placeholders
directory = cd; % default directory when prompting user to select a file

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
                VarToLoad = varargin{index+1};
                index = index + 2;
            case 'CMap'
                CMap = varargin{index+1};
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

if ~exist('saveFile', 'var') || isempty(Filename)
    [Filename, p] = uiputfile({'*.avi'},'Save file as',directory);
    if isnumeric(Filename)
        return
    end
    Filename = fullfile(p, Filename);
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
        temp = load(ImageFiles{findex}, VarToLoad, '-mat');
        Images{findex} = temp.(VarToLoad);
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
if ~isempty(StimID)
    showStimID = true;
elseif isempty(StimID) && showStimID
    StimID = 1:numStims;
end
if isnumeric(StimID)
    if isrow(StimID)
        StimID = StimID';
    end
    StimID = cellstr(num2str(StimID));
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
if ischar(CMap)
    switch CMap
        case 'HiLo'
            CMap = HiLoColormap(CLim(1), CLim(2));
        case 'b2r'
            CMap = b2r(CLim(1),CLim(2));
        case 'gray'
            CMap = gray(128);
        case 'red'
            CMap = [linspace(0,1,128)', zeros(128,1), zeros(128,1)];
        case 'green'
            CMap = [zeros(128,1), linspace(0,1,128)', zeros(128,1)];
        case 'blue'
            CMap = [zeros(128,1), zeros(128,1), linspace(0,1,128)'];
        case 'parula'
            CMap = parula(128);
    end
end

% Add stim marker to colormap
if showStimMarker && ~showStimID && ~showColorBar
    CMap = [CMap;[1,1,1]];
end

%% Save each stimulus average to video

% Open video
fprintf('Writing video: %s...', Filename);
vidObj = VideoWriter(Filename,'Motion JPEG AVI');
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
        
        if sindex == StimIndex(1) && findex == 1
            % Determine frame dimensions
            [H, W, ~] = size(Image);
            Dist = min(round(H/20),round(W/20));
        end
        
        % Display and save image to file
        if ~showStimID && ~showColorBar % don't need to display image as no overlays
            Image = (Image-CLim(1))/range(CLim); % scale by CLim
            if ~showStimMarker
                Image = round(Image*(size(CMap,1)-1)+1); % convert to colormap indexing
            else
                Image = round(Image*(size(CMap,1)-2)+1); % convert to colormap indexing
                if sindex ~= ControlIndex && findex>=StimFrameIndex(index,1) && findex<=StimFrameIndex(index,2)
                    Image(H-Dist*2:H-Dist,W-Dist*2:W-Dist) = size(CMap,1); % add stim marker
                end
            end
            Image = ind2rgb(Image, CMap);        % convert to RGB image
            writeVideo(vidObj, Image);           % save to file
            
        else % need to display image before saving
            
            if sindex == StimIndex(1) && findex == 1
                % Create figure
                hF = figure('Units', 'Pixels', 'Position', [50, 50, 1450, 950], 'Color', 'w');
                hA = axes('Parent', hF);
            end
            
            % Display Image
            imagesc(Image, CLim);
            set(hA,'DataAspectRatio',[1 1 1]);
            axis off; hold on;
            colormap(CMap);
            
            % Place Stimulus mark
            if showStimMarker && sindex ~= ControlIndex && findex>=StimFrameIndex(index,1) && findex<=StimFrameIndex(index,2)
                patch([W-Dist*2; W-Dist; W-Dist; W-Dist*2],...
                    [H-Dist*2; H-Dist*2; H-Dist; H-Dist],...
                    'white','EdgeColor','white');
            end
            
            % Place ID number
            if showStimID
                if sindex == ControlIndex
                    text(Dist, Dist, 'control', 'FontSize', 25, 'Color', 'w', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
                else
                    text(Dist, Dist, StimID{index}, 'FontSize', 25, 'Color', 'w', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
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