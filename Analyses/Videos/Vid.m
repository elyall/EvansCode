function SaveFile = Vid(SaveFile,Images,varargin)
%VID    Saves image frames to video file

% Default parameters that can be changed
Afilt = false;
% Afilt = fspecial('gaussian',5,1);
Tfilt = false;
% Tfilt = 1/100*ones(1,100);
frameRate = 15.45/4;
MCdata = {[]};
Maps = {[]};
MergeType = 'mean'; % 'mean' or 'blend'
StimIndex = false; % vector, logical or numeric
showColorBar = false;
Crop = false;
% Crop = [32.51, 0, 729.98, 512];
CLim = [];
CMapType = 'gray';

indMap = [];
Dim = [];
Map = [];

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'frameRate','FrameRate'}
                frameRate = varargin{index+1};
                index = index + 2;
            case 'StimIndex'
                StimIndex = varargin{index+1};
                index = index + 2;
            case 'showColorBar'
                showColorBar = true;
                index = index + 1;
            case 'Crop'
                Crop = varargin{index+1};
                index = index + 2;
            case 'CLim'
                CLim = varargin{index+1};
                index = index + 2;
            case 'Afilter'
                Afilt = varargin{index+1};
                index = index + 2;
            case 'Tfilter'
                Tfilt = varargin{index+1};
                index = index + 2;
            case 'CMapType'
                CMapType = varargin{index+1};
                index = index + 2;
            case 'MCdata'
                MCdata = varargin{index+1};
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
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Images','var') || isempty(Images)
    [Images,directory] = uigetfile({'*.sbx;*.tiff'}, 'Select image files:', directory,'MultiSelect','on');
    if isnumeric(Images)
        return
    end
    Images = fullfile(directory,Images);
end

if ~exist('SaveFile','var') || isempty(SaveFile)
    if iscellstr(Images)
        [~,SaveFile,~]=fileparts(Images{1});
    else
        SaveFile='vid';
    end
    [SaveFile,directory]=uiputfile({'*.avi'},'Save as',fullfile(directory,[SaveFile,'.avi']));
    if isempty(SaveFile)
        return
    end
    SaveFile=fullfile(directory,SaveFile);
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
        Images{findex} = load2P(Images{findex});
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
        Images{findex} = squeeze(Images{findex}(:,:,1,1,:));
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

% Convert maps to array
Maps = cat(1,Maps{:});

% Build composite image
if numFiles > 1
    if isempty(indMap)
        [Dim, refMap, indMap] = mapFoVs(Maps, 'type', MergeType);
    end
    temp = Images;
    Images = nan([refMap.ImageSize,numFrames]);
    for findex = 1:numFrames
        current = cellfun(@(x) x(:,:,findex), temp, 'UniformOutput', false);
        Images(:,:,findex) = createImage(current, Maps, 'speed', 'pretty', 'Map', indMap, Dim, refMap, 'OutputView', refMap);
    end
else
    Images = Images{1};
end


%% Determine color info

% Determine Color Limits
if isempty(CLim)
    temp = double(Images(:,:,round(end/2)));
    CLim = prctile(temp(:), [.01,99.99]);
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
    case 'hot'
        cmap = hot(128);
end


%% Determine stimuli info
if ~isequal(StimIndex,false) % Determine frame dimensions
    [H, W, ~] = size(Images);
    StimH = round(H/20);
    StimW = round(W/20);
end


%% Save each stimulus average to video
fprintf('Writing video: %s...', SaveFile);

% Open video
vidObj = VideoWriter(SaveFile,'Motion JPEG AVI');
set(vidObj, 'FrameRate', frameRate);
open(vidObj);

% Create figure
hF = figure('Units', 'Pixels', 'Position', [50, 50, 1450, 950], 'Color', 'w');
hA = axes('Parent', hF);

% Save each frame to video
parfor_progress(numFrames);
for findex = 1:numFrames
    
    % Display Image
    imagesc(Images(:,:,findex), CLim);
    axis equal off;
    colormap(cmap);
    
    % Place ID number
    if ~isequal(StimIndex,false) && StimIndex(findex)
        hold on;
        if islogical(StimIndex)
            patch([W-StimW*2; W-StimW; W-StimW; W-StimW*2],...
                [H-StimH*2; H-StimH*2; H-StimH; H-StimH],...
                'magenta','EdgeColor','magenta');
        else
            text(H-StimH, W-StimW, sprintf('%d',StimIndex(findex)), 'FontSize', 20, 'Color', 'm', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        end
    end
    
    % Display color bar
    if showColorBar
        cbH = colorbar;
        ylabel(cbH, 'Fluorescence (A.U.)');
    end
    
    % Write to video
    drawnow
    pause(0.01);
    frame = getframe(hF);
    writeVideo(vidObj, frame.cdata);
    
    parfor_progress;
end

close(vidObj);
close(hF);
parfor_progress(0);

fprintf('\tfinished\n');
