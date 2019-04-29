function Spikes = estimateSpikeTiming(Data)
% Data is numFrames x numROIs


% frameRate = 15.45;


%% Check input arguments
directory = cd;
if ~exist('Data','var') || isempty(Data)
    [Data, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(Data)
        return
    end
    Data = fullfile(p,Data);
end


%% Load ROI data
if ischar(Data)
    load(Data, 'ROIdata','-mat');
    Data = ROIdata;
end
if isstruct(Data)
    Data = gatherROIdata(Data,'rawdata')';
end
numROIs = size(Data,2);


%% Estimate spike timing

% % Initialize variables
% V.dt = 1/frameRate;  % time step size
% V.NCells = 1; % number of cells in each ROI

% Cycle through ROIs estimating spike timing
fprintf('Estimating spikes for %d neurons...', numROIs);
Spikes = zeros(size(Data,1),numROIs);
parfor rindex = 1:numROIs
    [~, Spikes(:,rindex), ~] = deconvolveCa(Data(:,rindex));
%     Spikes(:,rindex) = fast_oopsi(data, V);  
end
fprintf('\tComplete\n');

