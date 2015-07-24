function [SaveFile, Data, Map] = vidAverageWholeFoV(SaveFile, Data, IndicesToSave, Map, varargin)
% Data is output from this function, or cell array of strings of the
% ExperimentFiles to load data from.

% Display settings
type = 'AvgTrialdFoF';
CMapType = 'HiLo';
individualFiles = false;
showStimID = false;
showStimMarker = false;
showColorBar = false;
frameRate = 15.45;
speedUp = 0;
mergetype = 'quick'; % 'quick' or 'pretty'
Crop = false;
CLim = [];
blur = false;
h = fspecial('gaussian',5,1);

%% Check input arguments
if ~exist('SaveFile', 'var') || isempty(SaveFile)
    directory = CanalSettings('DataDirectory');
    [SaveFile, p] = uiputfile({'*.avi'},'Save file as',directory);
    if isnumeric(SaveFile)
        return
    end
    SaveFile = fullfile(p, SaveFile);
end

if ~exist('Data','var') || isempty(Data)
    directory = CanalSettings('ExperimentDirectory');
    [Data, p] = uigetfile({'*.mat'},'Choose Experiment file',directory,'MultiSelect','on');
    if isnumeric(Data)
        return
    end
    if iscellstr(Data)
        Data = fullfile(p, Data);
    else
        error('Use ''vidAverageStim'' instead for a single FoV');
    end
elseif ischar(Data)
    Data = {Data};
end

if ~exist('IndicesToSave', 'var') || isempty(IndicesToSave)
    IndicesToSave = 'all';
end

if ~exist('Map', 'var')
    Map = [];
end

index = 1;
while index<=length(varargin)
    switch varargin{index}
        case 'showStimID'
            showStimID = true;
            index = index + 1;
        case 'individualFiles'
            individualFiles = true;
            index = index + 1;
        case 'showStimMarker'
            showStimMarker = true;
            index = index + 1;
        case 'showColorBar'
            showColorBar = true;
            index = index + 1;
        case 'speedUp'
            speedUp = varargin{index+1};
            index = index + 2;
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
        case 'blur'
            blur = true;
            index = index + 1;
        case 'type'
            type = varargin{index+1};
            index = index + 2;
        case 'CMapType'
            CMapType = varargin{index+1};
            index = index + 2;
        otherwise
            warning('Argument ''%s'' not recognized',varargin{index});
            index = index + 1;
    end
end

%% Load In Data
numFiles = numel(Data);
fprintf('Will be saving average evoked dF/F to ''%s''\n', SaveFile);
if iscellstr(Data)
    fprintf(' \tloading %d files...\n', numFiles);
    ExperimentFiles = Data;
    clear Data;
    for findex = 1:numFiles
        Data(findex) = load(ExperimentFiles{findex}, 'AnalysisInfo', type, 'Map', '-mat');
        if ~isfield(Data, 'Map')
            Data(findex).Map = imref2d([size(Data(findex).AnalysisInfo, 1), size(Data(findex).AnalysisInfo, 2)]);
        end
        fprintf('\tloaded %d: %s\n', findex, ExperimentFiles{findex});
    end
    for findex = 1:numFiles
        Data(findex).filetype = ExperimentFiles{findex};        
        Data(findex).origMap = Data(findex).Map;
        Data(findex).cropped = false;
    end
end

ndim = ndims(Data(1).(type));
indexString = ':,:,';
for dindex = 3:ndim-1
    indexString = [indexString, '1,'];
end
indexString = [indexString, 'findex'];

%% Determine stimuli to save
numStims = numel(Data(1).(type));
if ischar(IndicesToSave) && strcmp(IndicesToSave, 'all')
    IndicesToSave = 1:numStims;
end

%% Determine stimulus periods
if showStimMarker
    stimIndex = zeros(numStims, 2);
    for sindex = 2:numStims % assumes first stim is the control trials
        stimIndex(sindex, :) = mode(Data(1).AnalysisInfo.TrialStimFrames(Data(1).AnalysisInfo.StimID==sindex-1,:), 1);
    end
end
    
%% Determine color limits
if isempty(CLim)
    findex = 30:50;
    temp = eval(sprintf('Data(1).(type){round(numStims/2)}(%s);', indexString));
    CLim = prctile(temp(:), [.01,99.99]);
end

%% Determine frame rate
if speedUp
    % frameRate = Data(1).frameRate*speedUp;
    frameRate = frameRate*speedUp;
end

%% Initialize Map for Generating Images
temp = cell(numFiles,1);
if isempty(Map) && numFiles > 1
    for findex = 1:numFiles
        temp{findex} = Data(findex).(type){4}(:,:,:,30);
    end
    [~,~,Map,Data] = createMultiFoVImage(Data, temp, mergetype, Crop);
    fprintf('Completed map\n');
%     assignin('base', 'Data', Data);
%     assignin('base', 'Map', Map);
end

%% Save each stimulus average to video

% Open video or determine filetypes to save to
if ~individualFiles
    fprintf('Writing video: %s\n', SaveFile);
    vidObj = VideoWriter(SaveFile,'Motion JPEG AVI');
    set(vidObj, 'FrameRate', frameRate);
    open(vidObj);
else
    [p,f,e] = fileparts(SaveFile);
    SaveFile = cell(numel(IndicesToSave), 1);
    for sindex = 1:numel(IndicesToSave)
        SaveFile{sindex} = fullfile(p,strcat(f,num2str(IndicesToSave(sindex)),e));
    end
end

hF = figure('Units', 'Pixels', 'Position', [50, 50, 1400, 800]);
hA = axes('Parent', hF);

switch CMapType
    case 'HiLo'
        cmap = HiLoColormap(CLim(1), CLim(2));
    case 'b2r'
        cmap = b2r(CLim(1),CLim(2));
    case 'gray'
        cmap = gray(128);
end

for sindex = IndicesToSave
    
    % For saving each stimulus to a file
    if individualFiles
        s = find(IndicesToSave == sindex);
        fprintf('\nWriting stim video: %s', SaveFile{s});
        vidObj = VideoWriter(SaveFile{s},'Motion JPEG AVI');
        set(vidObj, 'FrameRate', frameRate);
        open(vidObj);
    else
        fprintf('\n\twriting s%d:', sindex);
    end
    
    % Save each frame to file
    nFrames = size(Data(1).(type){1});
    nFrames = nFrames(end); %assumes frames are indexed last
    
    for findex = 1:nFrames
        
        % Build up frame
        for fileindex = 1:numFiles
            eval(sprintf('temp{fileindex} = Data(fileindex).(type){sindex}(%s);', indexString));
        end
        if numFiles > 1
            Image = createMultiFoVImage(Data, temp, mergetype, Crop, Map, 'Color', 'index');
        else
            Image = temp{1};
        end
        if sindex == IndicesToSave(1) && findex == 1
            [H, W, ~] = size(Image);
            StimH = round(H/20);
            StimW = round(W/20);
        end
        if blur
            Image = imfilter(Image, h);
        end
        
        % Display Image
        imagesc(Image, 'Parent', hA, CLim);
        axis off;
        colormap(cmap);
        
        % Place Stimulus mark
        if showStimMarker && findex>=stimIndex(sindex,1) && findex<=stimIndex(sindex,2)
            patch([W-StimW*2; W-StimW; W-StimW; W-StimW*2],...
                [H-StimH*2; H-StimH*2; H-StimH; H-StimH],...
                'white','EdgeColor','white');
        end
        
        % Place ID number
        if showStimID
            text(StimH, StimW, sprintf('%d', sindex), 'FontSize', 20, 'Color', 'cyan', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top');
        end
        
        % Display color bar
        if showColorBar
            cbH = colorbar;
            ylabel(cbH, 'dF/F');
        end
        
        % Write to video
        frame = getframe(hA);
        writeVideo(vidObj, frame.cdata);
        
        fprintf('\t%d', findex);
    end
    
    if individualFiles
        close(vidObj);
    end
end

if ~individualFiles
    close(vidObj);
end

close(hF);

fprintf('\nfinished\n');