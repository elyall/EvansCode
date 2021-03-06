function Images = loadWT(ImageFiles)

H = 512;
W = 640;

% H = 1024;
% W = 1280;

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
        ImageFiles = dir(ImageFiles);
        ImageFiles = fullfile(p, {ImageFiles(~[ImageFiles(:).isdir]).name});
    else
        ImageFiles = {ImageFiles};
    end
end

%% Initialize output
numFiles = numel(ImageFiles);
Images = zeros(H, W, numFiles);

%% Load in images
for findex = 1:numFiles
    
    % Load frame
    [~,~,ext] = fileparts(ImageFiles{findex});
    switch ext
        
        case '.raw'
            fid = fopen(ImageFiles{findex}, 'r');
            img = permute(reshape(fread(fid, inf, 'uint8=>uint8'), 3, W, H), [3,2,1]);
            Images(:,:,findex) = img(:,:,1);
            fclose(fid);
            
        case '.tif'
            img = imread(ImageFiles{findex}, 1);
            Images(:,:,findex) = img(:,:,1);
            
        case '.pgm'
            Images(:,:,findex) = imread(ImageFiles{findex});
            
    end
end