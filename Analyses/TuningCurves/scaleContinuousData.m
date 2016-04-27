function [Data, Bounds, numSamples] = scaleContinuousData(Data, Bounds, varargin)

numSamples = 100;
normalized = false;
FlipMap = false;
exponent = 1;

%% Parse input arguments
index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case 'numSamples'
                numSamples = varargin{index+1};
                index = index + 2;
            case 'exponent'
                exponent = varargin{index+1};
                index = index + 2;
            case 'normalized'
                normalized = true;
                index = index + 1;
            case 'flip'
                FlipMap = true;
                index = index + 1;
            otherwise
                warning('Argument ''%s'' not recognized',varargin{index});
                index = index + 1;
        end
    catch
        warning('Argument %d not recognized',index);
        index = index + 1;
    end
end

if ~exist('Bounds', 'var') || isempty(Bounds)
    Bounds = nan(1,2);
end


%% Determine Bounds
if isnan(Bounds(1))
    Bounds(1) = min(Data);
end
if isnan(Bounds(2))
    Bounds(2) = max(Data);
end

if Bounds(2) < Bounds(1)
    FlipMap = true;
    Bounds = flip(Bounds);
end

%% Scale data

% Cut off ends out of range
Data(Data < Bounds(1)) = Bounds(1);
Data(Data > Bounds(2)) = Bounds(2);

% Shift and scale range to between 0 to 1
Data = (Data - Bounds(1))/diff(Bounds);
if FlipMap
    Data = abs(Data-1);
end

% Scale data exponentially
Data = Data.^exponent;
Data = (Data - min(Data))/range(Data);

% Scale range so 1 now becomes the number of samples available
if ~normalized
    Data = round(Data*(numSamples-1)) + 1;
end

