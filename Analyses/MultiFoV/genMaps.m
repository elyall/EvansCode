function Maps = genMaps(Data)

Buffer = 1;
outputType = 'cell'; % 'array' or 'cell'

% Create maps
numFiles = numel(Data);
Dims = cellfun(@(x) [size(x,1),size(x,2)], Data, 'UniformOutput', false);
Dims = cat(1,Dims{:});
Maps = imref2d;
for f = 1:numFiles
    Maps(f) = imref2d(Dims(f,:));
end

% Determine panels dimensions
H = floor(sqrt(numFiles));  % number of panels in height
W = ceil(numFiles/H);       % number of panels in width
h = nan(H,1);
w = nan(W,1);
for index = 1:H
    h(index) = max(Dims((1:H:numFiles)+index-1,1)); % max height of current row
end
for index = 1:W
    w(index) = max(Dims(H*(index-1)+1:H*index,2)); % max width of current column
end
h = h+2*Buffer;
w = w+2*Buffer;

% Offset panels
for index = 1:numFiles
    [a,b] = ind2sub([H,W],index);
    Maps(index).YWorldLimits = Maps(index).YWorldLimits + sum(h(1:a-1)) + Buffer;
    Maps(index).XWorldLimits = Maps(index).XWorldLimits + sum(w(1:b-1)) + Buffer;
end

% Reformat output
if strcmp(outputType,'cell')
    temp = Maps;
    Maps = cell(numFiles,1);
    for findex = 1:numFiles
        Maps{findex} = temp(findex);
    end
end