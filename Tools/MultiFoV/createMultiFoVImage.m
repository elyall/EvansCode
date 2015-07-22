function [Image, Origin, Map, Data] = createMultiFoVImage(Data, Type, Speed, Crop, Map, varargin)
% Speed = 'quick' or 'pretty'
% Crop = false, true (for prompt), or numFiles x 4: xmin, ymin, width, height
% Map - only matters for 'pretty' mode

blur = false;
filt = fspecial('gaussian', 5, 1);
Color = 'redgreen'; % 'redgreen' or 'index'
Scaling = 'Independent'; % 'Independent' or 'Joint'

%% Parse input arguments
if ~exist('Data','var') || isempty(Data)
    directory = CanalSettings('DataDirectory');
    [Data, p] = uigetfile({'*.mat;*.ciexp'},'Choose Experiment file(s)', directory, 'MultiSelect', 'on');
    if isnumeric(Data)
        return
    end
    if iscell(Data)
        for findex = 1:numel(Data)
            Data{findex} = fullfile(p,Data{findex});
        end
    elseif ischar(Data)
        Data = {fullfile(p,Data)};
    end
elseif ischar(Data)
    Data = {Data};
end
numFiles = numel(Data);

if ~exist('Type', 'var') || isempty(Type)
    Type = 'average'; % 'average', 'variance', 'max', 'min', 'none'
end

if ~exist('Speed', 'var') || isempty(Speed)
    Speed = 'quick'; % 'quick' or 'pretty'
end

if ~exist('Crop', 'var') || isempty(Crop)
    Crop = false; % matrix (numFiles x 4:[xmin, ymin, width, height]) or 'true' for prompting
%     Crop = repmat([32.51, 0, 729.98, 512], numFiles, 1);
end
if ~islogical(Crop) && size(Crop, 1) == 1 && numFiles > 1
    Crop = repmat(Crop, numFiles, 1);
end

index = 1;
while index<=length(varargin)
    switch varargin{index}
        case 'blur'
            blur = true;
            index = index + 1;
        case 'filter'
            filter = varargin{index+1};
            index = index + 2;
        case 'Color'
            Color = varargin{index+1};
            index = index + 2;
        case 'Scaling'
            Scaling = varargin{index+1};
            index = index + 2;
        case 'files'
            ExperimentFiles = varargin{index+1};
            index = index + 2;
        otherwise
            warning('Argument ''%s'' not recognized',varargin{index});
            index = index + 1;
    end
end

%% Load in Data
if iscellstr(Data)
    ExperimentFiles = Data;
    clear Data;
    if ~iscell(Type) && any(strcmp(Type, {'average', 'variance', 'max', 'min'}))
        for findex = 1:numFiles
            Data(findex) = load(ExperimentFiles{findex}, 'Map', 'ImageFiles', '-mat');
            fprintf('\n\tloaded %d: %s', findex, ExperimentFiles{findex});
        end
    else
        for findex = 1:numFiles
            Data(findex) = load(ExperimentFiles{findex}, 'Map', '-mat');
        end
    end
    for findex = 1:numFiles
        Data(findex).file = ExperimentFiles{findex};
        Data(findex).origMap = Data(findex).Map;
        Data(findex).cropped = false;
    end
elseif isstruct(Data)
    if ~iscell(Type) && any(strcmp(Type, {'average', 'variance', 'max', 'min'})) && ~isfield(Data, 'ImageFiles')
        if isempty(ExperimentFiles)
            error('First argument is a ''Data'' struct, but does not contain necessary ''ImageFiles'' field. Please pass in the ''ExperimentFiles'' as an argument');
        else
            for findex = 1:numFiles
                load(ExperimentFiles{findex}, 'ImageFiles', '-mat');
                Data(findex).ImageFiles = ImageFiles;
                fprintf('\n\tloaded %d: %s', findex, ExperimentFiles{findex});
            end
        end
    end
end

%% Select and prepare images
for findex = 1:numFiles
    if iscell(Type)
        Data(findex).Image = Type{findex};
    else
        switch Type
            case 'average'
                Data(findex).Image = Data(findex).ImageFiles.Average;
            case 'variance'
                Data(findex).Image = Data(findex).ImageFiles.Var;
            case 'max'
                Data(findex).Image = Data(findex).ImageFiles.Max;
            case 'min'
                Data(findex).Image = Data(findex).ImageFiles.Min;
            case 'none'
%                 [H,W,~,~] = size(Data(findex).ImageFiles.Average);
                Data(findex).Image = zeros(512,796);
                Speed = 'quick';
                Color = 'index';
            case 'white'
%                 [H,W,~,~] = size(Data(findex).ImageFiles.Average);
                Data(findex).Image = ones(512,796);
                Speed = 'quick';
                Color = 'index';
        end
    end
    Data(findex).Image = permute(Data(findex).Image, [1,2,4,3]);
    Data(findex).Image = Data(findex).Image(:,:,:,1); % throw out all but the first z-slice
    [H, W, C] = size(Data(findex).Image);
    switch Color
        case 'redgreen'
            if C == 1
                Data(findex).Image = cat(3, zeros(H, W), Data(findex).Image, zeros(H, W)); % green only
                ChannelIndex = 2;
            elseif C == 2
                Data(findex).Image = cat(3, Data(findex).Image(:,:,[2,1]), zeros(H, W)); % green & red data
                ChannelIndex = [2,1];
            end
        case 'index'
            if C>1
                warning('Multiple channels present, only displaying first channel.');
                Data(findex).Image = Data(findex).Image(:,:,1);
            end
    end
    
    % Adjust Image
    Data(findex).current = Data(findex).Image;
    switch Color
        case 'redgreen'
            for cindex = 1:C
                % Determine color limits
                temp = Data(findex).Image(:,:,ChannelIndex(cindex));
                Data(findex).clim(cindex).limits = [min(temp(:)), max(temp(:))];
                Data(findex).clim(cindex).current = [Data(findex).clim(cindex).limits(1), prctile(temp(:),99.98)];
                % Scale colormap
                temp = (temp-Data(findex).clim(cindex).current(1))./(Data(findex).clim(cindex).current(2)-Data(findex).clim(cindex).current(1));
                temp(temp<0) = 0;
                temp(temp>1) = 1;
                Data(findex).current(:,:,ChannelIndex(cindex)) = temp;
            end
    end
    
    % Make sure each has a map
    if isempty(Data(findex).Map)
        Data(findex).Map = imref2d(size(Data(findex).current(:,:,1)));
    end
    
    % Crop image
    if (islogical(Crop) && Crop == true) || ~islogical(Crop)
        if islogical(Crop) && Crop == true
            [~, rect] = imcrop(Data(findex).current(:,:,min(2,end)));
        elseif ~islogical(Crop)
            [~, rect] = imcrop(Data(findex).current, Crop(findex, :));
        end
        rect([1,2]) = floor(rect([1,2]))+1;
        rect([3,4]) = ceil(rect([3,4]));
        Data(findex).current = Data(findex).current(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1,:);
        Data(findex).Map.XWorldLimits(1) = Data(findex).origMap.XWorldLimits(1)+rect(1);
        Data(findex).Map.XWorldLimits(2) = Data(findex).Map.XWorldLimits(1)+rect(3);
        Data(findex).Map.YWorldLimits(1) = Data(findex).origMap.YWorldLimits(1)+rect(2)-1;
        Data(findex).Map.YWorldLimits(2) = Data(findex).Map.YWorldLimits(1)+rect(4);
        Data(findex).Map.ImageSize = rect([4,3]);
        Data(findex).cropped = true;
    end

end

% Close cropping figure
if islogical(Crop) && Crop == true
    close gcf;
end

%% Build image
switch Color
    case 'redgreen'
        numChannels = 3;
    case 'index'
        numChannels = 1;
end

switch Speed
    
    case 'pretty'
        if ~exist('Map', 'var') || isempty(Map)
            [Data, Origin, Map] = buildCompositeMap(Data);
        else
            Origin = false;
        end
        [H, W, ~] = size(Map);
        Image = zeros(H, W, numChannels);
        for findex = 1:findex
            % Add average frame to composite image
            Image(Data(findex).ylim(1):Data(findex).ylim(2),Data(findex).xlim(1):Data(findex).xlim(2),:)...
                = Image(Data(findex).ylim(1):Data(findex).ylim(2),Data(findex).xlim(1):Data(findex).xlim(2),:)...
                + repmat(Map(Data(findex).ylim(1):Data(findex).ylim(2),Data(findex).xlim(1):Data(findex).xlim(2),findex),1,1,numChannels).*Data(findex).current;
        end
        
    case 'quick'
        for cindex = 1:size(Data(1).current, 3);
            temp = Data(1).current(:,:,cindex);
            Map = Data(1).Map;
            for findex = 2:numFiles
                [temp, Map] = imfuse(...
                    temp,...
                    Map,...
                    Data(findex).current(:,:,cindex),...
                    Data(findex).Map,...
                    'blend',...
                    'Scaling', Scaling);
            end
            switch Color
                case 'redgreen'
                    if cindex == 1
                        Image = cat(3, temp, zeros(size(temp,1),size(temp,2),2));
                    else
                        Image(:,:,cindex) = temp;
                    end
                case 'index'
                    Image = temp;
            end
            if ischar(Type) && strcmp(Type, 'white')
                Image(:) = 1;
            end
        end
        Origin = [Map.XWorldLimits(1), Map.YWorldLimits(1)];
end