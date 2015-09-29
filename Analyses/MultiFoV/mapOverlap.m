function Map = mapOverlap(MapsA, MapsB, varargin)
% Used to determine the region in which two datasets overlap

Crop = false;


%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'crop'
                Crop = varargin{index+1};
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

if ~exist('MapsA', 'var') || isempty(MapsA)
    [MapsA,p] = uigetfile({'*.exp;*.align'}, 'Select map files for first dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(MapsA)
        return
    elseif iscellstr(MapsA)
        MapsA = fullfile(p, MapsA);
    else
        MapsA = {fullfile(p, MapsA)};
    end
    directory = p;
elseif ischar(MapsA)
    MapsA = {MapsA};
end

if ~exist('MapsB', 'var') || isempty(MapsB)
    [MapsB,p] = uigetfile({'*.exp;*.align'}, 'Select map files for second dataset:', directory, 'MultiSelect', 'on');
    if isnumeric(MapsB)
        return
    elseif iscellstr(MapsB)
        MapsB = fullfile(p, MapsB);
    else
        MapsB = {fullfile(p, MapsB)};
    end
elseif ischar(MapsB)
    MapsB = {MapsB};
end


%% Load in maps from the two datasets
if iscellstr(MapsA)
    MapFiles = MapsA;
    MapsA = imref2d();
    for findex = 1:numel(MapFiles)
        load(MapFiles{findex}, 'Map', '-mat');
        MapsA(findex) = Map;
        clear Map;
    end
end

if iscellstr(MapsB)
    MapFiles = MapsB;
    MapsB = imref2d();
    for findex = 1:numel(MapFiles)
        load(MapFiles{findex}, 'Map', '-mat');
        MapsB(findex) = Map;
        clear Map;
    end
end

% Crop maps
if ~islogical(Crop) || Crop ~= false
    [~, MapsA] = crop([], Crop, MapsA);
    [~, MapsB] = crop([], Crop, MapsB);
end


%% Determine region of overlap

% Determine maps for each dataset
Maps = imref2d();
if numel(MapsA)>1
    [~, Maps(1)] = mapFoVs(MapsA, 'type', 'index');
else
    Maps(1) = MapsA;
end
if numel(MapsB)>1
    [~, Maps(2)] = mapFoVs(MapsB, 'type', 'index');
else
    Maps(2) = MapsB;
end

% Determine region of overlap between the two datasets
XLim = [-inf, inf];
YLim = [-inf, inf];
for findex = 1:2
    XLim(1) = max(XLim(1), Maps(findex).XWorldLimits(1));
    XLim(2) = min(XLim(2), Maps(findex).XWorldLimits(2));
    YLim(1) = max(YLim(1), Maps(findex).YWorldLimits(1));
    YLim(2) = min(YLim(2), Maps(findex).YWorldLimits(2));
end
H = diff(YLim);
W = diff(XLim);

% Create reference map for overlapping region
Map = imref2d([H,W], XLim, YLim);


