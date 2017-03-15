function [Image, Map, Images, Maps, indMap] = createImage(Images, Maps, varargin)


speed = 'quick'; % 'quick' or 'pretty'

% type = 'blend'; % 'mean' or 'blend' (for pretty only)
scaling = 'Independent'; % 'Independent' or 'Joint' (for quick only)

Crop = false;
% Crop = true; % matrix (numFiles x 4:[xmin, ymin, width, height]) or 'true' for prompting
% Crop = repmat([32.51, 0, 729.98, 512], 1, 1);

filt = false;
% filt = fspecial('gaussian', 5, 1);

filler = 'mean'; % 'mean' or scalar

outMap = []; % map of output view

imageType = 'average'; % for loading images only
directory = cd;

indMap = []; % can input blended map if previously created
Dim = [];
Map = [];

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'speed'
                speed = varargin{index+1};
                index = index + 2;
            case 'type'
                imageType = varargin{index+1};
                index = index + 2;
            case 'scaling'
                scaling = varargin{index+1};
                index = index + 2;
            case 'crop'
                Crop = varargin{index+1};
                index = index + 2;
            case 'filter'
                filt = varargin{index+1};
                index = index + 2;
            case 'OutputView'
                outMap = varargin{index+1};
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
    [Images, p] = uigetfile({'*.exp;*.mat'},'Choose Images files', directory, 'MultiSelect', 'on');
    if isnumeric(Images)
        return
    end
    Images = fullfile(p,Images);
    directory = p;
end
if ischar(Images)
    Images = {Images};
end

if ~exist('Maps','var') || isempty(Maps)
    [Maps, p] = uigetfile({'*.exp;*.mat'},'Choose map files', directory, 'MultiSelect', 'on');
    if isnumeric(Maps)
        return
    end
    Maps = fullfile(p,Maps);
end
if ischar(Maps)
    Maps = {Maps};
end

numFiles = numel(Images);


%% Load in maps
if iscellstr(Maps)
    MapFiles = Maps;
    Maps = imref2d();
    for findex = 1:numFiles
        load(MapFiles{findex}, 'Map', '-mat');
        Maps(findex) = Map;
    end
end


%% Load in images
if iscellstr(Images)
    ImageFiles = Images;
    for findex = 1:numFiles
        [~,~,ext] = fileparts(ImageFiles{findex});
        switch ext
            case '.exp'
                load(ImageFiles{findex}, 'ImageFiles', '-mat');
                switch imageType
                    case 'average'
                        Images{findex} = squeeze(ImageFiles.Average(:,:,1,1,:));
                end
            case '.align'
                load(ImageFiles{findex},'m','-mat');
                Images{findex} = m;
        end
    end
end

% Crop images
if ~islogical(Crop) || Crop ~= false
    [Images, Maps] = crop(Images, Crop, Maps);
end

% Filter images
if ~islogical(filt)
    for findex = 1:numFiles
        Images{findex} = imfilter(Images{findex}, filt);
    end
end

% Determine Image dimensions
if isrow(Images)
    Images = Images';
end
Dimensions = cell2mat(cellfun(@size, Images, 'UniformOutput', false));
if size(Dimensions,2) == 2
    Dimensions = cat(2, Dimensions, ones(numFiles,1));
end

% Scale images
% Max = cellfun(@(x) max(reshape(x,size(x,1)*size(x,2),size(x,3))), Images, 'UniformOutput', false);
% Max = max(cell2mat(Max));
% Min = cellfun(@(x) min(reshape(x,size(x,1)*size(x,2),size(x,3))), Images, 'UniformOutput', false);
% Min = min(cell2mat(Min));
% for findex = 1:numFiles
%     for cindex = 1:Dimensions(findex,3)
%         temp = (Images{findex}(:,:,cindex)-Min(cindex))./(Max(cindex)-Min(cindex));
%         temp(temp<0) = 0;
%         temp(temp>1) = 1;
%         Images{findex}(:,:,cindex) = temp;
%     end
% end


%% Build image
if numFiles > 1
    switch speed
        
        case {'pretty','blend'} % Build Pretty Image
            
            % Create map
            if isempty(indMap)
                [Dim, Map, indMap] = mapFoVs(Maps, 'type', 'blend');
            end
            [H,W,~] = size(indMap);
            
            % Create Image
            Image = zeros(H, W, max(Dimensions(:,3)));
            for findex = 1:numFiles
                Image(Dim(findex,2)+1:sum(Dim(findex,[2,4])),Dim(findex,1)+1:sum(Dim(findex,[1,3])),1:Dimensions(findex,3))...
                    = Image(Dim(findex,2)+1:sum(Dim(findex,[2,4])),Dim(findex,1)+1:sum(Dim(findex,[1,3])),1:Dimensions(findex,3))...
                    + repmat(indMap(Dim(findex,2)+1:sum(Dim(findex,[2,4])),Dim(findex,1)+1:sum(Dim(findex,[1,3])),findex),[1,1,Dimensions(findex,3)]).*Images{findex};
            end
            
            
        case {'quick','mean'} % Build Quick Image
            %             for cindex = 1:Dimensions(1,3);
            %                 Image = Images{1}(:,:,cindex);
            Image = Images{1}(:,:,1);
            Map = Maps(1);
            for findex = 2:numFiles
                [Image, Map] = imfuse(...
                    Image,...
                    Map,...
                    Images{findex}(:,:,1),...
                    Maps(findex),...
                    'blend',...
                    'Scaling', scaling);
            end
            Image = double(Image);
            %             end
            %         Origin = [Map.XWorldLimits(1), Map.YWorldLimits(1)];
            
            [~, ~, indMap] = mapFoVs(Maps, 'type', 'index'); % for determining what holes to fill
    end
    
    temp = Image(any(indMap,3));
    if ischar(filler)
        Image(~any(indMap,3)) = nanmean(temp(:));
        Image(isnan(Image)) = nanmean(temp(:));
    else
        Image(~any(indMap,3)) = filler;
        Image(isnan(Image)) = filler;
    end
    
else
    
    Image = Images{1};
    Map = Maps(1);
    
end


%% Shift image
if ~isempty(outMap)
    [Image, Map] = imwarp(Image, Map, affine2d(), 'OutputView', outMap);
end

