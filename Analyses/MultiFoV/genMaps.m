function Maps = genMaps(Data)

Buffer = [1,1];
H = 0;
W = 0;
outputType = 'array'; % 'array' or 'cell'


%% Initialize maps
if iscell(Data)
    numFiles = numel(Data);
    Dims = cellfun(@(x) [size(x,1),size(x,2)], Data, 'UniformOutput', false);
    Dims = cat(1,Dims{:});
elseif isnumeric(Data)
    N = ndims(Data);
    numFiles = size(Data,N);
    Dims = repmat([size(Data,1),size(Data,2)],numFiles,1);
end
Maps = imref2d;
for f = 1:numFiles
    Maps(f) = imref2d(Dims(f,:));
end


%% Determine panels dimensions
if ~H && ~W
    H = floor(sqrt(numFiles));  % number of panels in height
    W = ceil(numFiles/H);       % number of panels in width
elseif ~H
    H = ceil(numFiles/W);       % number of panels in height
elseif ~W
    W = ceil(numFiles/H);       % number of panels in width
end

% % Different dimension for each row/column
% Dims = cat(1,Dims,nan(H*W-numFiles,2));                  % add nans for reshaping (make square)
% h = nanmean(reshape(Dims(:,1), [W,H]), 1) + 2*Buffer(1); % height of each row
% w = nanmean(reshape(Dims(:,2), [W,H]), 2) + 2*Buffer(2); % width of each column

% Same height and width for each column and row respectively
h = repmat(max(Dims(:,1)) + 2*Buffer(1), [1,H]); % height of each row
w = repmat(max(Dims(:,2)) + 2*Buffer(2), [W,1]); % width of each column


%% Offset panels
[Y,X] = meshgrid(1:H,1:W);
Y = bsxfun(@times, Y-1, h)+Buffer(1); % calculate y offset for each panel
X = bsxfun(@times, X-1, w)+Buffer(2); % calculate x offset for each panel
% Index = reshape(1:H*W,W,H);           % transposed location of each panel
for ind = 1:numFiles
    [a,b] = ind2sub([W,H],ind);                               % determine current panel's location
%     [a,b] = find(Index==ind,1);                               % determine current panel's location
    Maps(ind).YWorldLimits = Maps(ind).YWorldLimits + Y(a,b); % add y offset for that panel
    Maps(ind).XWorldLimits = Maps(ind).XWorldLimits + X(a,b); % add x offset for that panel
end

% Format output
if strcmp(outputType,'cell')
    temp = Maps;
    Maps = cell(numFiles,1);
    for findex = 1:numFiles
        Maps{findex} = temp(findex);
    end
end