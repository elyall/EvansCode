function [Filename, CLim, indMap, Dim, Map] = vid(Filename, Images, varargin)
%vid Saves images to a video file
%   FILENAME = vid() prompts user to select a FILENAME to save to and then
%   select one or more image files to load data from. Each frame will be
%   merged across datasets, specified by the MAPS and MERGETYPE properties.
%
%   FILENAME = vidStim(FILENAME, IMAGES) saves IMAGES to the video file
%   FILENAME. IMAGES can be a cell array of filenames to load data from, or
%   a cell array 3-dimensional image data [H x W x F]. Each frame will be
%   merged across datasets. To concatenate datasets in time, input a cell
%   array of filenames as an element of the input cell array.

% TO DO: use insertText to display StimIDs without displaying to figure.


% Default parameters that can be adjusted

% Overlays & colorbar
StimFrameIndex = [];        % frame indices or logical vector of length numFrames specifying which frames to display a stimulus marker on
Text = {};                  % cell array of text to overlay
TextIndex = [];             % array of length numFrames indexing which text to display on each frame
Color = [1,1,1];            % color of overlays
showColorBar = false;       % whether to display the colorbar

% Specify data
MCdata = {[]};              % cell array of MCdata structures
Crop = false;               % false or [numFiles x 4] vector specifying number of pixels to remove from edges (see: crop)
borderLims = false;         % false or [numFiles x 4] vector specifying number of pixels along edges to set to zero, inclusive (top, bottom, left, right) (faster than crop)
Afilt = false;               % 2D filter to apply to each frame
% Afilt = fspecial('gaussian',5,1);               % 2D filter to apply to each frame
Tfilt = false;              % 1D filter to apply across time
% Tfilt = 1/100*ones(1,100);              % 1D filter to apply across time

% Merging properties
Maps = [];                  % cell array or array of imref2d objects specifying location of each input dataset
MergeType = 'mean';         % 'mean' or 'blend' specifying how to merge the datasets ('blend' takes a long time)
% Only used if want to skip repeating long blend process when making
% multiple videos (see: mapFoVs)
indMap = [];                % H x W x numFiles indexing array (see: mapFoVs)
Dim = [];                   % numFiles x 4 specifying the locations and size of each dataset (see: mapFoVs)
Map = [];                   % imref2d object of the composite data (see: mapFoVs)

% Display properties
CMap = 'parula';            % Nx3 colormap, or string specifying the colormap of the video
CLim = [];                  % 1x2 vector specifying the color limits ([] uses max and min)
frameRate = 15.45*2;        % scalar specifying the frame rate of the output video
pixelSize = [1,1];          % 1x2 vector specifying the aspect ratio of the Height and Width
outputSize = [];            % 1x2 vector specifying the desired Height and Width of the output video in pixels


% Placeholders
directory = cd; % default directory when prompting user to select a file

%% Check input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'StimFrameIndex'
                StimFrameIndex = varargin{index+1};
                index = index + 2;
            case 'Text'
                Text = varargin{index+1};
                index = index + 2;
            case 'TextIndex'
                TextIndex = varargin{index+1};
                index = index + 2;
            case 'Color'
                Color = varargin{index+1};
                index = index + 2;
            case 'showColorBar'
                showColorBar = true;
                index = index + 1;
            case 'MCdata'
                MCdata = varargin{index+1};
                index = index + 2;
            case 'Crop'
                Crop = varargin{index+1};
                index = index + 2;
            case 'borderLims'
                borderLims = varargin{index+1};
                index = index + 2;
            case 'Afilter'
                Afilt = varargin{index+1};
                index = index + 2;
            case 'Tfilter'
                Tfilt = varargin{index+1};
                index = index + 2;
            case 'Maps'
                Maps = varargin{index+1};
                index = index + 2;
            case 'MergeType'
                MergeType = varargin{index+1};
                index = index + 2;
            case 'Map'
                indMap = varargin{index+1};
                Dim = varargin{index+2};
                Map = varargin{index+3};
                index = index + 4;
            case 'CMap'
                CMap = varargin{index+1};
                index = index + 2;
            case 'CLim'
                CLim = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
                index = index + 2;
            case 'pixelSize'
                pixelSize = varargin{index+1};
                index = index + 2;
            case 'outputSize'
                outputSize = varargin{index+1};
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

if ~exist('Images', 'var') || isempty(Images)
    [Images,p] = uigetfile({'*.exp;*.align'}, 'Select image files for first dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(Images)
        return
    end
    Images = fullfile(p, Images);
end


%% Adjust inputs for number of files
if ~iscell(Images)
    Images = {Images};
end
numFiles = numel(Images);

if ~iscell(MCdata)
    MCdata = {MCdata};
end
if numel(MCdata)==1 && numFiles>1
    MCdata = repmat(MCdata, numFiles, 1);
end

if ~iscell(Maps)
    Maps = {Maps};
end
if numel(Maps)==1 && numFiles>1
    Maps = repmat(Maps, numFiles, 1);
end

if ~isequal(Crop,false) && size(Crop,1)==1 && numFiles>1
    Crop = repmat(Crop,numFiles,1);
end

%% Create and adjust images
for findex = 1:numFiles
    
    % Load images
    if iscellstr(Images{findex}) || ischar(Images{findex})
        Images{findex} = load2P(Images{findex},'Frames',[1,inf]);
    elseif iscell(Images{findex})
        Images{findex} = cat(ndims,Images{findex}{:});
    end
    
    % Motion correct images
    if ~isempty(MCdata{findex})
        if ischar(MCdata{findex})
            temp = load(MCdata{findex},'MCdata','-mat');
            MCdata{findex} = temp.MCdata;
        end
        Images{findex} = applyMotionCorrection(Images{findex}, MCdata{findex});
    end
    
    % Remove unwanted dimensions and convert to double
    if ndims(Images{findex})>3
        Images{findex} = squeeze(Images{findex});
    end
    Images{findex} = double(Images{findex});
    
    % Load map
    if ischar(Maps{findex})
        temp = load(Maps{findex},'Map','-mat');
        Maps{findex} = temp;
    elseif isempty(Maps{findex})
        Maps{findex} = imref2d([size(Images{findex},1),size(Images{findex},2)]);
    end
    
    % Crop images
    if ~isequal(Crop,false)
        [Images{findex},Maps{findex}] = crop(Images{findex}, Crop(findex,:), Maps{findex});
    end
    
    % Set border to zero
    if any(borderLims)
        Images{findex}(1:borderLims(1),:,:) = 0;            % top
        Images{findex}(end-borderLims(2)+1:end,:,:) = 0;    % bottom
        Images{findex}(:,1:borderLims(3),:) = 0;            % left
        Images{findex}(:,end-borderLims(4)+1:end,:) = 0;    % right
    end

    % Filter images in space
    if ~isequal(Afilt,false)
        Images{findex} = imfilter(Images{findex},Afilt);
    end
    
    % Filter images in time
    if ~isequal(Tfilt,false)
        Images{findex} = filter(Tfilt,1,double(Images{findex}),[],3);
    end
    
end
numFrames = size(Images{1},3);

if iscell(Maps)
    Maps = cat(1,Maps{:}); % convert maps to array
end

% Build composite images
if numFiles > 1
    if isempty(indMap)
        [Dim, refMap, indMap] = mapFoVs(Maps, 'type', MergeType); % generate indexing array
    end
    temp = Images;
    Images = nan([refMap.ImageSize,numFrames]); % initialize output array
    for findex = 1:numFrames
        current = cellfun(@(x) x(:,:,findex), temp, 'UniformOutput', false); % pull out current frame from each dataset
        Images(:,:,findex) = createImage(current, Maps, 'speed', MergeType, 'Map', indMap, Dim, refMap, 'OutputView', refMap); % build composite frame
    end
else
    Images = Images{1};
end

% Scale for aspect ratio and desired output size
sz = size(Images).*pixelSize; % determine how to scale to fix aspect ratio
if ~isempty(outputSize)
    sz = sz*min(outputSize./sz); % determine conversion to reach desired size while maintaining aspect ratio
end
Images = imresize(Images, sz); % scale images

% Add blanks to short edge to make output desired size
if ~isempty(outputSize)
    sz = size(Images);
    dimind = find(sz(1:2)-outputSize);
    num = outputSize(dimind)-sz(dimind);
    if dimind==1
        Images = cat(dimind,zeros(floor(num/2),sz(2),numFrames),Images,zeros(ceil(num/2),sz(2),numFrames));
    else
        Images = cat(dimind,zeros(sz(1),floor(num/2),numFrames),Images,zeros(sz(1),ceil(num/2),numFrames));
    end
end


%% Determine color info
if isempty(CLim)
    CLim = prctile(Images(:), [.01,99.99]);
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

% Add overlay color to colormap
if any(StimFrameIndex) || any(TextIndex)
    Overlay = true;
    CMap = [CMap;Color];
else
    Overlay = false;
end

%% Determine overlay info

% Determine frame dimensions
if any(StimFrameIndex)
    [H, W, ~] = size(Images);
    Dist = min(round(H/20),round(W/20));
    N = size(CMap,1)-2;
else
    N = size(CMap,1)-1;
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
        Image = createImage(temp, Maps, 'speed', MergeType, 'filter', filt, 'Map', indMap, Dim, Map);
        
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
                    text(Dist, Dist, Text{index}, 'FontSize', 25, 'Color', 'w', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
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


%% Save each stimulus average to video
fprintf('Writing video: %s...', SaveFile);

% Open video
vidObj = VideoWriter(SaveFile,'Motion JPEG AVI');
set(vidObj, 'FrameRate', frameRate);
open(vidObj);

% Create figure
if showColorbar
    hF = figure('Units', 'Pixels', 'Position', [50, 50, 1450, 950], 'Color', 'w');
    hA = axes('Parent', hF);
end

% Save each frame to video
parfor_progress(numFrames);
for findex = 1:numFrames
    Image = Images(:,:,findex);
    Image = (Image-CLim(1))/range(CLim); % scale by CLim
    
    % Add Stim Marker
    if ~showStimMarker
        Image = round(Image*N+1); % convert to colormap indexing
    else
        Image = round(Image*(size(CMap,1)-2)+1); % convert to colormap indexing
        if sindex ~= ControlIndex && findex>=StimFrameIndex(index,1) && findex<=StimFrameIndex(index,2)
            Image(H-Dist*2:H-Dist,W-Dist*2:W-Dist) = size(CMap,1); % add stim marker
        end
    end
    Image = ind2rgb(Image, CMap);        % convert to RGB image
    
    % Display Image
    if showColorbar
        image(Image);
        set(hA,'DataAspectRatio',[1 1 1]);
        axis off; hold on;
        colormap(CMap);
        cbH = colorbar;
        ylabel(cbH, 'Fluorescence (A.U.)');
        drawnow;
        frame = getframe(hF);
        img = frame.cdata;
    end
   
    % Write to video
    drawnow
    pause(0.01);
    writeVideo(vidObj, img);
    
    parfor_progress;
end

if showColorbar
    close(hF);
end
close(vidObj);
parfor_progress(0);

fprintf('\tfinished\n');