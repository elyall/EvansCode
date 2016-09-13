function [Data, Bounds, numSamples, tickLoc, tickLabel] = scaleContinuousData(Data, Bounds, varargin)

numSamples = 100;
normalized = false;
FlipMap = false;
exponent = 1;
numTicks = 5;

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
            case 'numTicks'
                numTicks = varargin{index+1};
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

if ~exist('Bounds', 'var') || isempty(Bounds)
    Bounds = nan(1,2);
end

if ~isa(Data,'double')
    Data = double(Data);
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

% % Scale data exponentially
% Data = Data.^exponent;
% Data = (Data - min(Data))/range(Data);

% Scale range so 1 now becomes the number of samples available
if ~normalized
    Data = round(Data*(numSamples-1)) + 1;
end


%% Determine ticks

% Determine tick distance
Range = range(Bounds);
unroundedTickSize = Range/(numTicks-1);
x = ceil(log10(unroundedTickSize)-1);
roundedTickRange = ceil(unroundedTickSize / 10^x) * 10^x;

% Determine ticks to show
tickBounds = [roundedTickRange*round(Bounds(1)/roundedTickRange), roundedTickRange*round(Bounds(2)/roundedTickRange)];
tickLabel = tickBounds(1):roundedTickRange:tickBounds(2);

% Determine tick locations
vals = Bounds(1):range(Bounds)/(numSamples-1):Bounds(2);
tickLoc = interp1(vals,1:numSamples,tickLabel,'linear','extrap');

% Old code before recognizing tick locations don't have to be integers
% tickCount = numel(Ticks);
% numValuesPerTick = round(numSamples/tickCount-1)-1;
% scale = linspace(tickBounds(1),tickBounds(2),numValuesPerTick*(tickCount-1)+1);
% increment = scale(2)-scale(1);
% newBounds = Bounds + [1,-1].*mod([-1,1].*Bounds,increment);
% roundTargets = newBounds(1):increment:newBounds(2);
% 
% [~,idx] = min(bsxfun(@(x,y)abs(x-y),Data',roundTargets')); %index of closest
% Data = roundTargets(idx); %extract values

