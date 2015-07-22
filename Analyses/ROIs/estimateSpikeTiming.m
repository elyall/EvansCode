function ROIdata = estimateSpikeTiming(ROIdata, NeuropilWeight, varargin)

saveOut = false;
saveFile = '';
frameRate = 15.5;

%% Check input arguments
if ~exist('ROIdata','var') || isempty(ROIdata)
    directory = CanalSettings('DataDirectory');
    [ROIdata, p] = uigetfile({'*.mat'},'Choose ROI file',directory);
    if isnumeric(ROIdata)
        return
    end
    ROIdata = fullfile(p,ROIdata);
end

if ~exist('NeuropilWeight', 'var')
    NeuropilWeight = 0.65;
end

index = 1;
while index<=length(varargin)
    try
        switch varargin{index}
            case {'save', 'Save'}
                saveOut = true;
                index = index + 1;
            case 'SaveFile'
                saveFile = varargin{index+1};
                index = index + 2;
            case 'frameRate'
                frameRate = varargin{index+1};
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


%% Load ROI data
if ischar(ROIdata)
    ROIFile = ROIdata;
    load(ROIdata, 'ROIdata');
    if saveOut && isempty(saveFile)
        saveFile = ROIFile;
    end
end
if isempty(saveFile)
    saveOut = false;
end
numROIs = numel(ROIdata.rois);

if ~isfield(ROIdata.rois, 'rawneuropil')
    NeuropilWeight = false;
end


%% Estimate spike timing

% Initialize variables
V.dt = 1/frameRate;  % time step size
V.NCells = 1; % number of cells in each ROI

% Cycle through ROIs estimating spike timing
fprintf('Estimating spikes for %d neurons...\n\tFinished ROI:', numROIs);
Spikes = zeros(numROIs, numel(ROIdata.rois(1).rawdata));
for rindex = 1:numROIs
    if NeuropilWeight
        Spikes(rindex,:) = fast_oopsi(ROIdata.rois(rindex).rawdata - NeuropilWeight*ROIdata.rois(rindex).rawneuropil, V);
    else
        Spikes(rindex,:) = fast_oopsi(ROIdata.rois(rindex).rawdata, V);
    end
    ROIdata.rois(rindex).rawspikes = Spikes(rindex,:);
    
    if ~mod(rindex, 5)
        fprintf('\t%d', rindex)
        if ~mod(rindex, 100)
            fprintf('\n\tFinished ROI:');
        end
    end
end
fprintf('\n');

%% Save to file
if saveOut
    save(saveFile, 'ROIdata', 'Spikes', '-append');
    fprintf('ROIdata saved to: %s\n', saveFile);
end
