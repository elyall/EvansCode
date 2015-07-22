function [Data, Origin, Map] = buildCompositeMap(varargin)
% input Maps
% output: H x W x N, where each element is the weight to be applied to the
% corresponding pixel in the Nth image being added to the composite

%% Parse input arguments
if numel(varargin) == 1 && isstruct(varargin{1})
    Data = varargin{1};
else
    Data = struct();
    for D = 1:numel(varargin)
        Data(D).Map = varargin{D};
    end
end
N = numel(Data);

%% Determine composite image size
XLim = [inf, -inf];
YLim = [inf, -inf];
for D = 1:N
    XLim(1) = min(XLim(1), Data(D).Map.XWorldLimits(1));
    XLim(2) = max(XLim(2), Data(D).Map.XWorldLimits(2));
    YLim(1) = min(YLim(1), Data(D).Map.YWorldLimits(1));
    YLim(2) = max(YLim(2), Data(D).Map.YWorldLimits(2));
end
H = diff(YLim);
W = diff(XLim);
Origin = [XLim(1),YLim(1)];

%% Determine location of images within map
Map = zeros(H, W, numel(Data));
center = [H/2, W/2];
for D = 1:N
    Data(D).ylim = [ceil(Data(D).Map.YWorldLimits(1)),floor(Data(D).Map.YWorldLimits(2))] - floor(YLim(1));
    Data(D).xlim = [ceil(Data(D).Map.XWorldLimits(1)),floor(Data(D).Map.XWorldLimits(2))] - floor(XLim(1));
    %     Data(D).center = [Data(D).ylim(1)+diff(Data(D).ylim)/2, Data(D).xlim(1)+diff(Data(D).xlim)/2];
    [~,Data(D).yCornerIndex] = min(abs(Data(D).ylim-center(1)));
    [~,Data(D).xCornerIndex] = min(abs(Data(D).xlim-center(1)));
    Data(D).corner = [Data(D).ylim(Data(D).yCornerIndex), Data(D).xlim(Data(D).xCornerIndex)];
    Map(Data(D).ylim(1):Data(D).ylim(2),Data(D).xlim(1):Data(D).xlim(2),D) = 1;
end

%% Average overlaying points
if nargout > 2
    fullMap = sum(Map, 3);
    for D = 1:N
        % Find data points current data set is overlapping with others
        [I,J] = find(fullMap(Data(D).ylim(1):Data(D).ylim(2),Data(D).xlim(1):Data(D).xlim(2)) > 1);
        J = J + Data(D).xlim(1)-1;
        I = I + Data(D).ylim(1)-1;
        % Determine image's borders
        left = [(Data(D).ylim(1):Data(D).ylim(2))', repmat(Data(D).xlim(1),Data(D).Map.ImageExtentInWorldY,1)];
        right = [(Data(D).ylim(1):Data(D).ylim(2))', repmat(Data(D).xlim(2),Data(D).Map.ImageExtentInWorldY,1)];
        top = [repmat(Data(D).ylim(2),Data(D).Map.ImageExtentInWorldX,1),(Data(D).xlim(1):Data(D).xlim(2))'];
        bottom = [repmat(Data(D).ylim(1),Data(D).Map.ImageExtentInWorldX,1),(Data(D).xlim(1):Data(D).xlim(2))'];
        borders = [left;right;top;bottom];
        % Determine distance of overlapping pixels from border
        dist = ipdm([I,J], borders, 'Subset','NearestNeighbor', 'Result', 'Structure');
        % Record distance from border as pixel's weight
        Map(sub2ind(size(Map),I,J,repmat(D,numel(I),1))) = [dist(:).distance];
    end
    % Normalize weights across overlapping images to sum to 1 (data closer
    % to the border will count less than data further from the border)
    Map = Map./repmat(sum(Map,3),1,1,N);
end

