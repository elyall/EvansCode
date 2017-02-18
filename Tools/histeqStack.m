function Images = histeqStack(Images,varargin)

NumTiles = [16,16];
Distribution = 'Exponential';

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'NumTiles'
                NumTiles = varargin{index+1};
                index = index + 2;
            case 'Distribution'
                Distribution = varargin{index+1};
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
    [Images,p] = uigetfile({'*.tif;*.Sbx'}, 'Select files containing maps:', directory);
    if isnumeric(Maps)
        return
    end
    Images = fullfile(p,Images);
end

%% Perform adaptive histogram equalization
if ischar(Images)
    Images = load2P(Images, 'Frames', [2,inf]);
end
while ndims(Images) < 5
    Images = permute(Images, [1:ndims(Images)-1, ndims(Images)+1, ndims(Images)]);
end

[~,~,D,C,F] = size(Images);
for d = 1:D
    for c = 1:C
        for f = 1:F
            temp = Images(:,:,d,c,f);
            Images(:,:,d,c,f) = adapthisteq(temp/max(temp(:)),'NumTiles',NumTiles,'Distribution',Distribution);
        end
    end
end



