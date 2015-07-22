function [Data, numSamples, CLim] = scaleContinuousData(Data, numSamples, lowerBound, upperBound)

if ~exist('lowerBound', 'var') || isempty(lowerBound)
    lowerBound = min(Data);
end

if ~exist('upperBound', 'var') || isempty(upperBound)
    upperBound = max(Data);
end

% Cut off ends out of range
Data(Data > upperBound) = upperBound;
Data(Data < lowerBound) = lowerBound;

% Shift and scale range to between 0 to 1
Data = (Data - lowerBound)/(upperBound-lowerBound);
Data = round(Data*numSamples);

% Fix issue with very bottom of range
Data(Data == 0) = 1;

% Define CLim
CLim = [lowerBound, upperBound];