function framecounter = loadFrameCounter(ImageFiles)

offset = 4;
directory = 'C:\Users\LightFieldPC\Documents\WhiskerTracking';

%% Parse input arguments
if ~exist('ImageFiles', 'var') || isempty(ImageFiles)
    [ImageFiles, p] = uigetfile({'*.raw;*.tif;*.pgm'}, 'Select Whisker Tracking Image Files', directory, 'MultiSelect', 'on');
    if isnumeric(ImageFiles)
        return
    elseif iscell(ImageFiles)
        for findex = 1:numel(ImageFiles)
            ImageFiles{findex} = fullfile(p, ImageFiles{findex});
        end
    elseif ischar(ImageFiles)
        ImageFiles = {fullfile(p, ImageFiles)};
    end
elseif ischar(ImageFiles)
    if isdir(ImageFiles) % select all files in input directory
        p = ImageFiles;
        ImageFiles = dir(p);
        ImageFiles = fullfile(p, {ImageFiles(~[ImageFiles(:).isdir]).name});
    else
        ImageFiles = {ImageFiles};
    end
end

%% Initialize output
numFiles = numel(ImageFiles);
framecounter = nan(numFiles, 4);

%% Load in images
for findex = 1:numFiles
    
    % Open file
    fid = fopen(ImageFiles{findex}, 'r');
    fseek(fid, offset, 'bof');
    
    % Read framecounter info
    framecounter(findex,:) = fread(fid, 4, 'ubit8');
    
    % Close file
    fclose(fid);
    
end