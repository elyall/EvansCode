function timestamps = loadTimeStamps(ImageFiles)
%(f,1)= second (0 to 127)
%(f,2)= cycle_count (0 to 7999)
%(f,3)= cycle_offset (0 to x -> x depends on configuration)

offset = 0;
directory = 'C:\Users\LightFieldPC\Documents\WhiskerTracking';

%% Parse input arguments
if ~exist('ImageFiles', 'var') || isempty(ImageFiles)
    [ImageFiles, p] = uigetfile({'*.raw;*.tif;*.pgm;*.avi'}, 'Select Whisker Tracking Image Files', directory, 'MultiSelect', 'on');
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
timestamps = nan(numFiles, 3);

%% Load in images
for findex = 1:numFiles
    
    % Open file
    fid = fopen(ImageFiles{findex}, 'r');
    fseek(fid, offset, 'bof');
    
    % Read timestamp info
    timestamps(findex, 1) = fread(fid, 1, 'ubit7');     % second (0 to 127)
    timestamps(findex, 2) = fread(fid, 1, 'ubit13');    % cycle_count (0 to 7999)
    timestamps(findex, 3) = fread(fid, 1, 'ubit12');    % cycle_offset (0 to x -> x depends on configuration)
    
%     bin = fread(fid, 32, 'ubit1');
%     timestamps(findex, 1) = binaryVectorToDecimal(bin(1:7)');
%     timestamps(findex, 2) = binaryVectorToDecimal(bin(8:20)');
%     timestamps(findex, 3) = binaryVectorToDecimal(bin(21:32)');
    
    % Close file
    fclose(fid);
    
end