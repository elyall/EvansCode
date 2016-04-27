function [Dim, refMap, indMap, Maps] = mapFoVs(Maps, varargin)
% Dim - first two values are X,Y distance from top left corner. second
% two values are width and height.

type = 'index'; % 'index' or 'mean' or 'blend'(indMap only)

directory = cd;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'Type','type'}
                type = varargin{index+1};
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

if ~exist('Maps', 'var')
    [Maps,p] = uigetfile({'*.exp;*.align'}, 'Select files containing maps:', directory, 'MultiSelect', 'on');
    if isnumeric(Maps)
        return
    elseif iscellstr(Maps)
        Maps = fullfile(p, Maps);
    else
        Maps = {fullfile(p, Maps)};
    end
elseif ischar(Maps)
    Maps = {Maps};
end

numFiles = numel(Maps);


%% Load in maps
if iscellstr(Maps)
    MapFiles = Maps;
    Maps = imref2d();
    for findex = 1:numFiles
        temp = load(MapFiles{findex}, 'Map', '-mat');
        if isfield(temp,'Map')
            Maps(findex) = temp.Map;
        end
        clear temp;
    end
end


%% Determine location of images
XLim = [inf, -inf];
YLim = [inf, -inf];
for findex = 1:numFiles
    XLim(1) = min(XLim(1), Maps(findex).XWorldLimits(1));
    XLim(2) = max(XLim(2), Maps(findex).XWorldLimits(2));
    YLim(1) = min(YLim(1), Maps(findex).YWorldLimits(1));
    YLim(2) = max(YLim(2), Maps(findex).YWorldLimits(2));
end
H = diff(YLim);
W = diff(XLim);


%% Determine image dimensions
Dim = zeros(numFiles, 4);
for findex = 1:numFiles
    Dim(findex, :) = [Maps(findex).XWorldLimits(1) - XLim(1), Maps(findex).YWorldLimits(1) - YLim(1), flip(Maps(findex).ImageSize)];
end


%% Reference map
if nargout > 1
    refMap = imref2d([H,W], XLim, YLim);
end


if nargout > 2
    %% Index map
    indMap = zeros(H, W, numFiles);
    for findex = 1:numFiles
        ylim = [ceil(Maps(findex).YWorldLimits(1)),floor(Maps(findex).YWorldLimits(2))] - floor(YLim(1));
        xlim = [ceil(Maps(findex).XWorldLimits(1)),floor(Maps(findex).XWorldLimits(2))] - floor(XLim(1));
        indMap(ylim(1):ylim(2),xlim(1):xlim(2),findex) = true;
    end
    
    
    %% Blend map
    switch type
        
        case 'index'
            indMap = logical(indMap);
            
        case 'mean'
            indMap = bsxfun(@rdivide, indMap, sum(indMap,3));
            
        case 'blend'
            fullMap = sum(indMap, 3);
            for findex = 1:numFiles
                
                % Find data points current data set is overlapping with others
                currentRegion = reshape(fullMap(logical(indMap(:,:,findex))), Dim(findex,4), Dim(findex,3));
                [I,J] = find(currentRegion > 1);
                J = J + Dim(findex,1);
                I = I + Dim(findex,2);
                
                % Determine image's borders
                lft = [(Dim(findex,2)+0.5:sum(Dim(findex,[2,4]))+.5)',          repmat(Dim(findex,1)+.5,Dim(findex,4)+1,1)];
                rgt = [(Dim(findex,2)+0.5:sum(Dim(findex,[2,4]))+.5)',          repmat(sum(Dim(findex,[1,3]))+.5,Dim(findex,4)+1,1)];
                top = [repmat(Dim(findex,2)+.5,Dim(findex,3)+1,1),              (Dim(findex,1)+.5:sum(Dim(findex,[1,3]))+.5)'];
                btm = [repmat(sum(Dim(findex,[2,4]))+.5,Dim(findex,3)+1,1),     (Dim(findex,1)+.5:sum(Dim(findex,[1,3]))+.5)'];
                borders = [lft;rgt;top;btm];
                
                % SANITY CHECK
                %subplot(1,numFiles,findex); imagesc(fullMap); hold on;
                %plot(lft(:,2),lft(:,1), 'r', 'LineWidth', 2);
                %plot(rgt(:,2),rgt(:,1), 'r', 'LineWidth', 2);
                %plot(top(:,2),top(:,1), 'r', 'LineWidth', 2);
                %plot(btm(:,2),btm(:,1), 'r', 'LineWidth', 2);
                
                % Determine distance of overlapping pixels from border
                dist = ipdm([I,J], borders, 'Subset','NearestNeighbor', 'Result', 'Structure');
                
                % Record distance from border as overlapping pixel's weight
                indMap(sub2ind([H,W,numFiles],I,J,repmat(findex,numel(I),1))) = [dist(:).distance];
                
            end
            
            % Normalize weights across overlapping regions to sum to 1
            % (data closer to the border will count less than data further
            % from the border)
            indMap = bsxfun(@rdivide, indMap, sum(indMap,3));
    end
end





