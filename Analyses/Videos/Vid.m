function saveFile = Vid(saveFile,Images,varargin)

% Filters
Afilt = false;
% Afilt = fspecial('gaussian',5,1);
Tfilt = false;
% Tfilt = 1/100*ones(1,100);

% Speed
frameRate = 15.45;

% Motion correction
MCdata = [];

% Display info
StimIndex = false; % vector, logical or numeric
showColorBar = false;
% Crop = false;
Crop = [32.51, 0, 729.98, 512];

% Color info
CLim = [];
CMapType = 'green';

% Loading only
Channel = 1;
Frames = [1 inf];

directory = cd;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'Frames'
                Frames = varargin{index+1};
                index = index + 2;
            case 'Channel'
                Channel = varargin{index+1};
                index = index + 2;
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
    [Images,directory] = uigetfile({'*.sbx;*.tiff'}, 'Select image files:', directory);
    if isnumeric(Images)
        return
    end
    Images = fullfile(directory,Images);
end

if ~exist('saveFile','var') || isempty(saveFile)
    if iscellstr(Images)
        [~,saveFile,~]=fileparts(Images{1});
    else
        saveFile='vid';
    end
    [saveFile,directory]=uiputfile({'*.avi'},'Save as',fullfile(directory,[saveFile,'.avi']));
    if isempty(saveFile)
        return
    end
    saveFile=fullfile(directory,saveFile);
end


%% Load Images
if iscellstr(Images) || ischar(Images)
    [Images, loadObj] = load2P(Images,'Frames',Frames,'Channel',Channel);
end
numFrames = size(Images,5);

% Motion correct images
if ~isempty(MCdata)
    Images = applyMotionCorrection(Images, MCdata, loadObj);
end

% Crop images
if ~isequal(Crop,false)
    Images = crop(Images, Crop);
end

% Filter images in space
if ~isequal(Afilt,false)
    Images = imfilter(Images,Afilt);
end
        
% Filter images in time
if ~isequal(Tfilt,false)
    Images = filter(Tfilt,1,double(Images),[],3);
end


%% Determine color info

% Determine Color Limits
if isempty(CLim)
    temp = double(Images(:,:,:,:,round(end/2)));
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
end


%% Determine stimuli info
if ~isequal(StimIndex,false) % Determine frame dimensions
    [H, W, ~] = size(Images);
    StimH = round(H/20);
    StimW = round(W/20);
end


%% Save each stimulus average to video
fprintf('Writing video: %s...', saveFile);

% Open video
vidObj = VideoWriter(saveFile,'Motion JPEG AVI');
set(vidObj, 'FrameRate', frameRate);
open(vidObj);

% Create figure
hF = figure('Units', 'Pixels', 'Position', [50, 50, 1450, 950], 'Color', 'w');
hA = axes('Parent', hF);

% Save each frame to video
parfor_progress(numFrames);
for findex = 1:numFrames
    
    % Display Image
    imagesc(Images(:,:,:,1,findex), CLim);
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
